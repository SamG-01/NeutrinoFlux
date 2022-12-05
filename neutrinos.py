from NeutrinoFlux.__init__ import *
from NeutrinoFlux.cross_sections import CrossSection, cross_section_GR, GR_bounds
from NeutrinoFlux.differential_flux import astro_flux, atmo_flux

class Neutrino:
    """Class containing the properties of a neutrino, including its cross sections and differential flux laws."""

    def __init__(self, flavor, anti) -> None:
        """Defines the basic properties of a particular neutrino. Of note, the cross sections for the neutrino can be modified before being passed into the event_rate function in main.py."""

        # Basic Properties
        self.flavor = flavor # choices: e, tau, mu
        self.anti = anti # choices: True, False

        # Name: for accessing atmo flux function
        self.string = "_" + ("anti" if anti else "") + "nu" + flavor

        # Cross Section Functions: dictionary. can be modified
        self.cross_sections = {
            "nc": cross_sections[anti]["nc"],
            "cc": cross_sections[anti]["cc"]
        }

        if anti and self.flavor == "e":
            self.cross_sections["GR"] = cross_sections["GR"]

        # Differential flux functions (default)
        self.flux_funcs = {
            "atmo": atmo_flux,
            "astro": astro_flux
        }

    def effective_volume(self, E, theta):
        """Returns the effective volume of the detector for a given neutrino at a certain energy and zenith angle."""
        return (1e3)**3 # m^3; 1 km^3
    
    def attenuation(self, E, theta):
        """Returns the absorption coefficient at a given energy and angle for a given cross section."""
        cs = self.cross_sections
        
        exponent = -1 * attenuation_parameter(theta)
        exponent *= sum(sigma(E) * sigma.targets_per_gram_earth for sigma in cs.values())
        return np.exp(exponent)
        
    def diff_flux(self, E, theta, flux_type, kwargs):
        """Returns the differential flux Phi(E, theta) for a given neutrino and flux type (astro or atmo).
        
        If astro, there are additional arguments for the spectral index gamma and the normalization phi_astro.
        
        If atmo, there are additional arguments for the month (January or July) and the flux_source (total, pi, k, pr, or conv)."""

        kwargs["self"] = self # for when the flux depends on the neutrino type

        if flux_type == "total":
            return sum(flux(E, theta, kwargs) for flux in self.flux_funcs.values())
        else:
            flux_func = self.flux_funcs[flux_type]
            return flux_func(E, theta, kwargs)


# Some default neutrino objects for testing purposes
default_neutrinos = {
    "nu_e": Neutrino("e", False),
    "nubar_e": Neutrino("e", True),
    "nu_mu": Neutrino("mu", False),
    "nubar_mu": Neutrino("mu", True),
    "nu_tau": Neutrino("tau", False),
    "nubar_tau": Neutrino("tau", True)
}
