from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np
from scipy.interpolate import interpn

from MCEq.core import MCEqRun
import mceq_config
import crflux.models as pm

from . import constants as C

_atmo_flux_path = Path(__file__).parent/"data/atmo_flux"


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
                    month: str = "January") -> np.ndarray:

        # MCEq defines theta = 0 when neutrinos are coming down from the atmosphere,
        # and theta > 90 when they are up-going through the earth.
        # We have the opposite coordinate system, so we need to convert our angle first.

        theta = 180 - theta
        cls.mceq_run.set_density_model(('MSIS00_IC', ('SouthPole', month)))

        # suppresses annoying logs from mceq set_theta
        with redirect_stdout(NullIO()):
            cls.mceq_run.set_theta_deg(180 - theta)
        cls.mceq_run.solve()
        return cls.mceq_run.get_solution(nu.atmo_name(source))

    @classmethod
    @C.ureg.wraps(C.flux_units, [None, "deg", None, None, None])
    def solve_flux(cls, theta: float,
                   nu: C.NeutrinoData,
                   source: str = "total",
                   month: str | list[str] = "Average") -> np.ndarray:
        """Solves for neutrino flux on the E grid for theta."""

        if month == "Average":
            month = C.months

        # averages the flux over the months given
        if not isinstance(month, str):
            return sum(cls.solve_flux(theta, nu, source, m)
                       for m in month)/len(month)

        assert source in cls.sources
        atmo_name = nu.atmo_name(source)

        # checks the class cache first
        key = (atmo_name, month)
        flux = cls.flux_cache.get(key, None)
        if flux is not None:
            return flux

        # if not, checks the data files
        data_file = _atmo_flux_path/month.lower()/atmo_name
        if data_file.exists():
            flux = np.loadtxt(str(data_file))
        else:
            # otherwise, solve for it directly
            flux = cls._solve_flux(theta, nu, source, month)

            # saves the data to a file
            # for future reference
            np.savetxt(str(data_file), flux)

        # adds the file to the class cache
        atmo_flux.flux_cache[key] = flux
        return flux

    @classmethod
    @C.ureg.wraps(C.flux_units, [None, "GeV", "deg", None, None, None, None])
    def interpolate(cls, E: float, theta: float,
                    nu: C.NeutrinoData,
                    source: str = "total",
                    month: str | list[str] = "January",
                    method: str = "cubic") -> float:
        """Interpolates the atmospheric flux for theta and E."""

        # solves for flux on the grid
        flux_data = cls.solve_flux(cls.theta_grid, nu, source, month)

        # expands any sparse grids
        E = E * np.ones_like(theta)
        theta = theta * np.ones_like(E)

        xi = np.stack((np.ravel(theta), np.ravel(E)), axis=-1)
        flux = interpn((cls.theta_grid.m, cls.E_grid.m),
                       flux_data.m, xi, method=method)

        # if input data is 0D, return the element
        if flux.size == 1:
            return flux.item()
        return flux.reshape(np.shape(E))
