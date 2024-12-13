from dataclasses import dataclass
from pint import Quantity, UnitRegistry

ureg = UnitRegistry(auto_reduce_dimensions=True)
ureg.setup_matplotlib(True)

pi = 3.141592653589793
months = ['January', 'February', 'March', 'April', 'May', 'June',
          'July', 'August', 'September', 'October', 'November', 'December']

# earth properties
R_E = 6371 * ureg.km
rho_earth_avg = 5.51 * ureg.g/ureg.cm**3

# used to convert mass density to number density
nucleon_concentration_earth = (1/ureg.amu).to("1/gram")

# number of electrons = number of protons, which make up about half the nucleon mass
electron_concentration_earth = 0.5 * nucleon_concentration_earth

# water properties; there are 10 electrons per 18 nucleons
nucleon_concentration_water = nucleon_concentration_earth
electron_concentration_water = 10/18 * nucleon_concentration_water

rho_ice = 0.917 * ureg.g/ureg.cm**3

# glashow resonance constants
G_F = 1.16638e-5 * ureg.GeV**-2
M_W = 80.379 * ureg.GeV/ureg.c**2

G_W = 0.02634
R_W_mu = 0.1057

sigma_0 = (G_F * M_W * ureg.hbar * ureg.c**3)**2 / pi

# differential flux
flux_units = (ureg.GeV * ureg.cm**2 * ureg.s * ureg.sr)**-1

# Astro Flux
C_0 = 3e-18 * flux_units
E_0 = 100 * ureg.TeV

gamma = 2.53  # works for each neutrino flavor
phi_astro = 1.66  # combined for nu + nubar

flavors = ("e", "tau", "mu")
anti = (True, False)


@dataclass
class NeutrinoData:
    flavor: str = "e"
    anti: bool = False

    def __post_init__(self):
        assert self.flavor in ["e", "tau", "mu"]
        assert isinstance(self.anti, bool)

    def __str__(self) -> str:
        return "nu" + ("bar" if self.anti else "") + "_" + self.flavor

    def atmo_name(self, source: str) -> str:
        return source + "_" + str(self).replace("_", "").replace("nubar", "antinu")

    def crsc_name(self, cc_or_nc: str = "cc") -> str:
        return "nu" + ("bar" if self.anti else "") + "_" + cc_or_nc
