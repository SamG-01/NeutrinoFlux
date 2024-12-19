from functools import partial

import numpy as np

from . import constants as C
from .cross_sections import CrossSection
from .differential_flux import astro_flux, atmo_flux
from .earth_density import attenuation_parameter


class Neutrino(C.NeutrinoData):
    """Class containing the properties of a neutrino,
    including its cross sections and differential flux laws."""

    def __post_init__(self) -> None:

        super().__post_init__()

        self.cross_sections: dict[str, CrossSection] = {
            cc_or_nc: CrossSection.from_cache(self.crsc_name(cc_or_nc))
            for cc_or_nc in ["cc", "nc"]
        }

        if self.anti and self.flavor == "e":
            self.cross_sections["gr"] = CrossSection.from_cache("gr")

        self.flux_funcs = {
            "astro": partial(astro_flux, nu=self),
            "atmo": partial(atmo_flux.interpolate, nu=self)
        }

    def effective_volume(self, E: C.Quantity, theta: C.Quantity) -> C.Quantity:
        """Returns the effective volume of the detector for a
        given neutrino at a certain energy and zenith angle."""

        return 1 * C.ureg.km**3

    def crsc_hits(self, E: C.Quantity, target: str = "earth",
                  cross_sections: dict[str, CrossSection] | None = None) -> C.Quantity:
        """Sums sigma(E) * target_conc for each cross section."""

        if cross_sections is None:
            cross_sections = self.cross_sections

        return sum(sigma.times_conc(E, target) for sigma in cross_sections.values())

    def earth_attenuation(self, E: C.Quantity, theta: C.Quantity,
                          cross_sections: dict[str, CrossSection] | None = None) -> C.Quantity:
        """Returns the earth absorption coefficient at a given energy
        and angle for the given cross sections."""

        exponent = -attenuation_parameter(theta)
        exponent = exponent * self.crsc_hits(E, "earth", cross_sections)
        return np.exp(exponent)

    def diff_flux(self, E: C.Quantity, theta: C.Quantity, flux_type: str, **kwargs) -> C.Quantity:
        """Returns the differential flux Phi(E, theta) for a given neutrino and flux type (astro or atmo).

        If astro, there are additional arguments for the spectral
        index gamma and the normalization phi_astro.

        If atmo, there are additional arguments for the month
        and the flux_source (total, pi, k, pr, or conv)."""

        kwargs["nu"] = self

        if flux_type == "total":
            return sum(flux_func(E, theta, **kwargs)
                       for flux_func in self.flux_funcs.values())

        flux_func = self.flux_funcs[flux_type]
        return flux_func(E, theta, **kwargs)


# Some default neutrino objects for testing purposes
default_neutrinos = {
    "nu_e": Neutrino("e", False),
    "nubar_e": Neutrino("e", True),
    "nu_mu": Neutrino("mu", False),
    "nubar_mu": Neutrino("mu", True),
    "nu_tau": Neutrino("tau", False),
    "nubar_tau": Neutrino("tau", True)
}
