from .imports import *
from .cross_sections import cross_section_GR
from .differential_flux import astro_flux, atmo_flux

class Neutrino():
    def __init__(self, flavor, anti) -> None:
        self.flavor = flavor # choices: e, tau, mu
        self.anti = anti # choices: True, False

        # Cross Sections        
        self.sigma = cross_section_total[anti]

        if flavor == "e" and anti:
            self.GR = cross_section_GR
        else:
            self.GR = lambda E: 0

        # For astro flux
        self.string = "total_" + ("anti" if anti else "") + "nu" + flavor

    def dN_dE(self, E, theta, flux_type, month):
        if flux_type == "astro":
            return astro_flux(E)
        elif flux_type == "atmo":
            return atmo_flux(E, theta, self.string, month)
        elif flux_type == "total":
            return astro_flux(E) + atmo_flux(E, theta, self.string, month)
        else:
            raise TypeError
