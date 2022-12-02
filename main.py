from NeutrinoFlux.__init__ import *
from NeutrinoFlux.neutrinos import Neutrino
from NeutrinoFlux.cross_sections import GR_a, GR_b

def attenuation(E, theta, sigma):
    """Returns the absorption coefficient at a given energy and angle for a given cross section."""
    exponent = attenuation_parameter(theta)
    exponent *= -1 * sum(cs.eval(E) * cs.targets_per_gram_earth for cs in sigma.values())
    return np.exp(exponent)

def integrand(E, theta, diff_flux, flux_type, sigma, atten, kwargs):
    """The integrand sigma * Phi * the attenuation factor."""
    integrand = mass_density_ice * sum(cs.eval(E) * cs.targets_per_gram_water for cs in sigma.values())
    integrand *= diff_flux(E, theta, flux_type, kwargs) * np.sin(theta)

    if atten:
        integrand *= attenuation(E, theta, sigma)
    
    return integrand

def event_rate(nu, flux_type, E_bounds, theta_bounds, diff_flux_kwargs=None, atten=True, GR_only=False):
    """Returns the yearly rate of neutrinos in the detector for a given neutrino type and flux type. If working with atmo flux, also can be configured based on the month of year and the source of the atmospheric flux."""
    
    if diff_flux_kwargs == None: diff_flux_kwargs = {}

    if GR_only:
        sigma = {"GR": nu.sigma["GR"]}
    else:
        sigma = nu.sigma

    args = (nu.diff_flux, flux_type, sigma, atten, diff_flux_kwargs)
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
    yearly_rate = 2 * np.pi * V_eff * rate_seconds * seconds_per_year

    return yearly_rate
