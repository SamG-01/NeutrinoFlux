from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
from scipy.interpolate import interpn

from MCEq.core import MCEqRun
import mceq_config
import crflux.models as pm

from . import constants as C

path = Path(__file__).parent/"data/atmo_flux"

class NullIO(StringIO):
    def write(self, txt):
        pass

# Sets energy bounds (in GeV)
mceq_config.e_min = 1e4
mceq_config.e_max = 1e12


def astro_flux(E: C.Quantity, *args,
               gamma: float = C.gamma,
               phi_astro: float = C.phi_astro,
               **kwargs) -> C.Quantity:
    """Flux law for astro neutrinos."""

    # astro flux for nu and nubar is approximately the same, so we divide phi_astro by 2
    return (C.C_0 * phi_astro / 2 * (E/C.E_0)**(-gamma)).to(C.flux_units)


class atmo_flux:
    """Flux law for atmospheric neutrinos."""

    sources = ["total", "pi", "k", "pr", "conv"]

    mceq_run = MCEqRun(
        # provide the string of the interaction model
        interaction_model='SIBYLL2.3c',

        # primary cosmic ray flux model

        # support a tuple (primary model class (not instance!), arguments)
        primary_model=(pm.HillasGaisser2012, "H3a"),

        # density model
        density_model=('MSIS00_IC', ('SouthPole', 'January')),

        # Zenith angle in degrees. 0=vertical, 90=horizontal
        theta_deg=0,
    )

    # Define equidistant grid in cos(theta)
    theta_grid = np.arccos(np.linspace(1, -1, 81) * C.ureg.km**0).to("deg")
    E_grid = mceq_run.e_grid * C.ureg.GeV

    flux_cache = {}

    @classmethod
    @np.vectorize(excluded=[0, 2, 3, 4], signature="()->(n)")
    def _solve_flux(cls, theta: float,
                    nu: C.NeutrinoData,
                    source: str = "total",
                    month: str | list[str] = "January"):

        if month == "Average":
            month = C.months

        if isinstance(month, str):
            assert source in cls.sources
        else:
            return np.mean([cls._solve_flux(theta, nu, source, m)
                           for m in month], axis=0)

        # supresses annoying logging from mceq set_theta
        with redirect_stdout(NullIO()):
            cls.mceq_run.set_theta_deg(theta)
        cls.mceq_run.solve()

        return cls.mceq_run.get_solution(nu.atmo_name(source))

    @classmethod
    @C.ureg.wraps(None, [None, "deg", None, None, None])
    def solve_flux(cls, theta: float,
                   nu: C.NeutrinoData,
                   source: str = "total",
                   month: str | list[str] = "January"):
        """Solves for neutrino flux on the E grid for theta."""

        if month == "Average":
            month = C.months

        if isinstance(month, str):
            assert source in cls.sources
        else:
            return sum(cls.solve_flux(theta, nu, source, m) for m in month)/len(month)

        cls.mceq_run.set_density_model(('MSIS00_IC', ('SouthPole', month)))
        return cls._solve_flux(theta, nu, source, month)

    @classmethod
    @C.ureg.wraps(C.flux_units, [None, "GeV", "deg", None, None, None, None])
    def interpolate(cls, E: float, theta: float,
                    nu: C.NeutrinoData,
                    source: str = "total",
                    month: str | list[str] = "January",
                    method: str = "cubic") -> float:
        """Interpolates the atmospheric flux for theta and E."""

        # MCEq defines theta = 0 when neutrinos are coming down from the atmosphere,
        # and theta > 90 when they are up-going through the earth.
        # We have the opposite coordinate system, so we need to convert our angle first.

        # solves for flux on the grid
        key = (nu.atmo_name(source), month)
        flux_data = cls.flux_cache.get(key, None)
        if flux_data is None:
            flux_data = cls.solve_flux(cls.theta_grid, nu, source, month)
            cls.flux_cache[key] = flux_data

        E = E * np.ones_like(theta)
        theta = theta * np.ones_like(E)

        xi = np.stack((180 - np.ravel(theta), np.ravel(E)), axis=-1)
        flux = interpn((cls.theta_grid.m, cls.E_grid.m),
                       flux_data, xi, method=method)

        # if input data is 0D, return the element
        if flux.size == 1:
            return flux.item()
        return flux.reshape(E.shape)

neutrinos = [C.NeutrinoData(flavor, anti) for flavor in C.flavors for anti in C.anti]
for m in C.months:
    mpath = path/m.lower()
    for neutrino in neutrinos:
        for flux_source in atmo_flux.sources:
            file = mpath/neutrino.atmo_name(flux_source)
            if not file.exists():
                flux_ = atmo_flux.solve_flux(atmo_flux.theta_grid, neutrino, flux_source, m)
                np.savetxt(str(file), flux_)
            else:
                atmo_flux.flux_cache[file.name, m] = np.loadtxt(str(file))
