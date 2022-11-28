from NeutrinoFlux.imports import *
from NeutrinoFlux.neutrinos import Neutrino
from NeutrinoFlux.cross_sections import GR_a, GR_b

def integrand(E, theta, diff_flux, flux_type, sigma, sigma_GR, attenuation, month, atmo_source):
    """The integrand sigma * dN/dE * the attenuation factor."""
    integrand = sigma(E) + 10/18 * sigma_GR(E)
    integrand *= diff_flux(E, theta, flux_type, month, atmo_source) * np.sin(theta)

    if attenuation:
        exponent = attenuation_function(theta)
        exponent *= -1 * (sigma(E) + sigma_GR(E)/2)
        integrand *= np.exp(exponent)
    
    return integrand

def event_rate(nu, flux_type, E_bounds, theta_bounds, atmo_source="total", month="January", attenuation=True, GR_only=False):
    """Returns the yearly rate of neutrinos in the detector for a given neutrino type and flux type. If working with atmo flux, also can be configured based on the month of year and the source of the atmospheric flux."""
    
    if GR_only: nu.sigma = lambda E: 0

    args = (nu.diff_flux, flux_type, nu.sigma, nu.sigma_GR, attenuation, month, atmo_source)
    E_a, E_b = E_bounds

    # Integrates the GR range of 4-8 PeV separately
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

    # Integrates over phi, and multiplies by other conversion factors to get the yearly rate.
    yearly_rate = 2 * np.pi * number_density_water * V_eff * rate_seconds * seconds_per_year

    return yearly_rate
