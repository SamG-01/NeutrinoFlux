from .imports import *

# Astro Flux
C0 = 3e-18 * 1/C.giga * 1/C.centi**2 # 1/(eV * m**2 * second * steradian)
E0 = 100 * C.tera #* eV
gamma = ufloat(2.53, 0.07) # works for each neutrino flavor
phi_astro = ufloat(1.66, 0.27) # combined for nu + nubar

def astro_flux(E):
    """Power law for astro neutrino flux."""
    return C0 * phi_astro.n / 2 * (E/E0)**(-gamma.n)

# Atmospheric Flux


def atmo_flux(E, theta, name):
    """Power law for atmospheric neutrino flux."""
    pass
