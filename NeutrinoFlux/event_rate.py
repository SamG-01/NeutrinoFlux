import numpy as np

from . import constants as C
from .cross_sections import CrossSection
from .neutrinos import Neutrino

def integrand(nu: Neutrino, E: C.Quantity, theta: C.Quantity,
              flux_type: str = "total", attenuate: bool = True,
              cross_sections: CrossSection | dict[str, CrossSection] | None = None,
              **flux_kwargs) -> C.Quantity:
    """The integrand rho * V_eff(E) * sigma(E) * Phi(E, theta) * attenuation(E, theta)."""

    ice_mass = C.rho_ice * nu.effective_volume(E, theta)

    if cross_sections is None:
        cross_sections = nu.cross_sections
    if isinstance(cross_sections, CrossSection):
        crsc_hits = cross_sections.times_conc(E, "ice")
    else:
        crsc_hits = nu.crsc_hits(E, "ice", cross_sections)

    Phi = nu.diff_flux(E, theta, flux_type, **flux_kwargs)
    attenuation = nu.earth_attenuation(E, theta) if attenuate else 1

    expr = ice_mass * crsc_hits * Phi * attenuation
    return expr

def event_rate(nu: Neutrino, E_bounds: C.Quantity | None = None,
               theta_bounds: C.Quantity | None = None,
               phi_bounds: C.Quantity | None = None,
               flux_type: str = "total",
               attenuate: bool = True,
               cross_sections: CrossSection | dict[str, CrossSection] | None = None,
               N: int = 1000, **flux_kwargs) -> C.Quantity:

    if E_bounds is None:
        E_bounds = [1e4, 1e12] * C.ureg.GeV
    if theta_bounds is None:
        theta_bounds = [0, np.pi] * C.ureg.rad
    if phi_bounds is None:
        phi_bounds = [0, 2*np.pi] * C.ureg.rad

    E_ = np.geomspace(*E_bounds.to("GeV").m, N) * C.ureg.GeV
    theta_ = np.linspace(*theta_bounds.to("rad"), N)
    E, theta = np.meshgrid(E_, theta_, sparse=True)

    integrand_ = integrand(nu, E, theta, flux_type, attenuate, cross_sections, **flux_kwargs)
    integral = np.trapz(np.trapz(integrand_ * np.sin(theta), theta_, axis=0), E_, axis=0)
    rate = (phi_bounds[1] - phi_bounds[0]) * integral

    return rate.to("1/year")
