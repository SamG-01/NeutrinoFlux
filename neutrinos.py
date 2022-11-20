from NeutrinoFlux.imports import *
from NeutrinoFlux.cross_sections import cross_section_GR
from NeutrinoFlux.differential_flux import astro_flux, atmo_flux

class Neutrino():
    def __init__(self, flavor, anti) -> None:
        # Basic Properties
        self.flavor = flavor # choices: e, tau, mu
        self.anti = anti # choices: True, False

        # Cross Sections   
        self.sigma_nc = cross_sections[anti]["nc"]
        self.sigma_cc = cross_sections[anti]["cc"]
        self.sigma = cross_sections[anti]["tot"]

        if flavor == "e" and anti:
            self.sigma_GR = cross_section_GR
        else:
            self.sigma_GR = lambda E: 0

        # For astro flux
        self.string = "_" + ("anti" if anti else "") + "nu" + flavor

    def dN_dE(self, E, theta, flux_type, month, atmo_source="total"):
        if flux_type == "astro":
            return astro_flux(E)
        elif flux_type == "atmo":
            return atmo_flux(E, theta, self.string, month, atmo_source)
        elif flux_type == "total":
            return astro_flux(E) + atmo_flux(E, theta, month, self.string, atmo_source)
        else:
            raise TypeError

default_neutrinos = [
    Neutrino("e", False), # nu_e
    Neutrino("e", True), # nubar_e
    Neutrino("mu", False), # nu_mu
    Neutrino("mu", True), # nubar_mu
    Neutrino("tau", False), # nu_tau
    Neutrino("tau", True) # nubar_tau
]
