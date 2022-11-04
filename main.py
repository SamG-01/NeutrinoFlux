from .imports import *
from .neutrinos import Neutrino
from .cross_sections import GR_a, GR_b

def integrand(E, theta, dN_dE, flux_type, sigma, GR, attenuation):
    """The integrand sigma * dN/dE * the attenuation factor."""
    integrand = sigma(E) + 10/18 * GR(E)
    integrand *= dN_dE(E, theta, flux_type) * np.sin(theta)

    if attenuation:
        #exponent = avg_rho_earth * x(theta)
        exponent = attenuation_function(theta)
        exponent *= -1 * (sigma(E) + GR(E))
        integrand *= np.exp(exponent)
    
    return integrand

def event_rate(flavor, anti, flux_type, E_bounds, theta_bounds, attenuation=True, GR_only=False):
    """Returns the yearly rate of neutrinos in the detector for a given neutrino type and flux type."""
    nu = Neutrino(flavor, anti)
    if GR_only: nu.sigma = lambda E: 0

    sigma, GR, dN_dE = nu.sigma, nu.GR, nu.dN_dE
    args = (dN_dE, flux_type, sigma, GR, attenuation)
    E_a, E_b = E_bounds

    if E_a <= GR_a and GR_b <= E_b:
        rate_seconds = nquad(
                    integrand, [(E_a, GR_a), theta_bounds], args=args
                )[0] + nquad(
                    integrand, [(GR_a, GR_b), theta_bounds], args=args
                )[0] + nquad(
                    integrand, [(GR_b, E_b), theta_bounds], args=args
                )[0]

    elif E_a > GR_a and GR_b <= E_b:
        rate_seconds = nquad(
                    integrand, [(E_a, GR_b), theta_bounds], args=args
                )[0] + nquad(
                    integrand, [(GR_b, E_b), theta_bounds], args=args
                )[0]

    elif E_a <= GR_a and E_b <= GR_b:
        rate_seconds = nquad(
                    integrand, [(E_a, GR_a), theta_bounds], args=args
                )[0] + nquad(
                    integrand, [(GR_a, E_b), theta_bounds], args=args
                )[0]

    else:
        rate_seconds = nquad(
                    integrand, [(E_a, E_b), theta_bounds], args=args
                )[0]

    yearly_rate = 2 * np.pi * number_density_water * V_eff * rate_seconds * seconds_per_year

    return yearly_rate
