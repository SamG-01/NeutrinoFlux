from NeutrinoFlux.__init__ import *
from NeutrinoFlux.cross_sections import GR_a, GR_b

def integrand(E, theta, nu, sigmas, flux_type, diff_flux_kwargs, atten):
    """The integrand V_eff(E) * sigma(E) * Phi(E, theta) * attenuation(E, theta)."""

    integrand = mass_density_ice * nu.effective_volume(E)
    integrand *= sum(sigma(E) * sigma.targets_per_gram_water for sigma in sigmas.values())
    integrand *= nu.diff_flux(E, theta, flux_type, diff_flux_kwargs)

    if atten: integrand *= nu.attenuation(E, theta)

    return integrand * np.sin(theta)

def event_rate(nu, flux_type, E_bounds, theta_bounds, phi_bounds=(0,2*np.pi), diff_flux_kwargs=None, atten=True, GR_only=False):
    """Returns the yearly rate of neutrinos in the detector for a given neutrino type and flux type. The entries in the diff_flux_kwargs dict are passed into the neutrino's diff_flux function and thus will vary based on the type of flux being computed."""
    
    if diff_flux_kwargs == None:
        diff_flux_kwargs = {}

    if GR_only:
        sigmas = {"GR": nu.cross_sections["GR"]}
    else:
        sigmas = nu.cross_sections

    args = (nu, sigmas, flux_type, diff_flux_kwargs, atten)
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
    phi_a, phi_b = phi_bounds
    yearly_rate = (phi_b - phi_a) * rate_seconds * seconds_per_year

    return yearly_rate
