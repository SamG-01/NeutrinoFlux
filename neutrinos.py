from NeutrinoFlux.imports import *
from NeutrinoFlux.cross_sections import CrossSection, cross_section_GR, GR_bounds
from NeutrinoFlux.differential_flux import astro_flux, atmo_flux

class Neutrino():
    """Class containing the properties of a high-energy neutrino, including its cross sections and differential flux laws."""
    def __init__(self, flavor, anti) -> None:
        """Defines the basic properties of a particular neutrino. Of note, the cross sections for the neutrino can be modified before being passed into the event_rate function in main.py."""

        # Basic Properties
        self.flavor = flavor # choices: e, tau, mu
        self.anti = anti # choices: True, False

        # Name: for accessing atmo flux function
        self.string = "_" + ("anti" if anti else "") + "nu" + flavor

        # Cross Section Functions: dictionary. can be modified
        self.sigma = {
            "nc": cross_sections[anti]["nc"],
            "cc": cross_sections[anti]["cc"]
        }
        if flavor == "e" and anti:
            self.sigma["GR"] = cross_sections["GR"]
        else:
            self.sigma["GR"] = CrossSection("blank", anti, lambda E: 0, True, 0, 0)


    def diff_flux(self, E, theta, flux_type, month, atmo_source="total"):
        """Returns the differential flux Phi(E, theta) for a given neutrino and flux type (astro or atmo). If atmo, arguments for the month (January or July) and the source (total, pi, k, pr, or conv)."""

        if flux_type == "astro":
            return astro_flux(E)
        elif flux_type == "atmo":
            return atmo_flux(E, theta, month, self.string, atmo_source)
        elif flux_type == "total":
            return astro_flux(E) + atmo_flux(E, theta, month, self.string, atmo_source)
        else:
            raise TypeError

# Some default neutrino objects for testing purposes
default_neutrinos = {
    "nu_e": Neutrino("e", False),
    "nubar_e": Neutrino("e", True),
    "nu_mu": Neutrino("mu", False),
    "nubar_mu": Neutrino("mu", True),
    "nu_tau": Neutrino("tau", False),
    "nubar_tau": Neutrino("tau", True)
}
