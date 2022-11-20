from NeutrinoFlux.imports import *
from NeutrinoFlux.neutrinos import default_neutrinos

flux_units = 1/C.giga * 1/C.centi**2 # 1/(eV * m**2 * second * steradian)

# Astro Flux
C0 = 3e-18 * flux_units
E0 = 100 * C.tera #* eV
gamma = ufloat(2.53, 0.07) # works for each neutrino flavor
phi_astro = ufloat(1.66, 0.27) # combined for nu + nubar

def astro_flux(E):
    """Power law for astro neutrino flux."""
    return C0 * phi_astro.n / 2 * (E/E0)**(-gamma.n)

# Atmospheric Flux
if __name__ == "__main__":
    mceq_run = MCEqRun(
        #provide the string of the interaction model
        interaction_model='SIBYLL2.3c',
        
        #primary cosmic ray flux model
        
        #support a tuple (primary model class (not instance!), arguments)
        primary_model = (pm.HillasGaisser2012, "H3a"),
        
        #density model
        density_model=('MSIS00_IC',('SouthPole','January')),

        # Zenith angle in degrees. 0=vertical, 90=horizontal
        theta_deg=0,
    )

    atmo_flux_sources = ["total", "pi", "k", "pr", "conv"]

    mceq_config.e_min = 1e13/C.giga
    mceq_config.e_max = 1e21/C.giga

    mceq_run.solve()

    #obtain energy grid (fixed) of the solution for the x-axis of the plots
    e_grid = mceq_run.e_grid

    #Define equidistant grid in cos(theta)
    angles = np.arccos(np.linspace(1,-1,20))*180./np.pi

    #Dictionary for results
    flux_funcs = {}

    for month in ("January", "July"):
        flux_funcs[month] = {}
        mceq_run.set_density_model(('MSIS00_IC',('SouthPole',month)))
        for nu in default_neutrinos:
            name = nu.string
            flux_funcs[month] = {name: {}}
            flux_data = {flux_source: [] for flux_source in atmo_flux_sources}

            for theta in angles:
                print(theta)
                mceq_run.set_theta_deg(theta)
                mceq_run.solve()
                for flux_source in atmo_flux_sources:
                    flux_data[flux_source].append(mceq_run.get_solution(flux_source + nu.string) * flux_units)

            for flux_source in atmo_flux_sources:
                flux_funcs[month][name][flux_source] = RectBivariateSpline(angles, e_grid, flux_data[flux_source])

    with open("Fits/atmo_flux", "wb") as f:
        pickle.dump(flux_funcs, f)

def atmo_flux(E, theta, month, neutrino_name, flux_source):
    """Power law for atmospheric neutrino flux."""
    theta = 180 * (1 - theta/np.pi)

    flux_func = flux_funcs[month][neutrino_name][flux_source]
    flux = flux_func(theta, E/C.giga)
    return flux[0][0]
