from NeutrinoFlux.imports import *
#from NeutrinoFlux.neutrinos import default_neutrinos

flux_units = 1/C.giga * 1/C.centi**2 # 1/(eV * m**2 * second * steradian)

# Astro Flux
C0 = 3e-18 * flux_units
E0 = 100 * C.tera #* eV
gamma = ufloat(2.53, 0.07) # works for each neutrino flavor
phi_astro = ufloat(1.66, 0.27) # combined for nu + nubar

def astro_flux(E):
    """Flux law for astro neutrinos."""
    # astro flux for nu and nubar is approximately the same, so we divide phi_astro by 2
    return C0 * phi_astro.n / 2 * (E/E0)**(-gamma.n)

# Atmospheric Flux
if __name__ == "__main__": # Performs fitting for atmo flux

    # First, sets up the interaction model
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

    # The different sources of atmospheric flux
    atmo_flux_sources = ["total", "pi", "k", "pr", "conv"]

    # Sets energy bounds
    mceq_config.e_min = 1e13/C.giga
    mceq_config.e_max = 1e21/C.giga

    mceq_run.solve()

    #obtain energy grid (fixed) of the solution for the x-axis of the plots
    e_grid = mceq_run.e_grid * C.giga

    #Define equidistant grid in cos(theta)
    angles = np.arccos(np.linspace(1,-1,20))*180./np.pi

    #Dictionary for results
    flux_funcs = {}

    for month in ("January", "July"):
        flux_funcs[month] = {}
        mceq_run.set_density_model(('MSIS00_IC',('SouthPole',month)))
        for anti in (True, False):
            for flavor in ("e", "mu", "tau"):
                name = "_" + ("anti" if anti else "") + "nu" + flavor
                flux_funcs[month][name] = {}
                flux_data = {flux_source: [] for flux_source in atmo_flux_sources}

                for theta in angles:
                    print(theta)

                    # Solves the system for a given theta
                    mceq_run.set_theta_deg(theta)
                    mceq_run.solve()

                    # Extracts solution for each type of neutrino, for each source of flux, at that theta
                    for flux_source in atmo_flux_sources:
                        flux_data[flux_source].append(mceq_run.get_solution(flux_source + name) * flux_units)

                # Interpolates the grid of solutions in theta, E into full functions for each flux source
                for flux_source in atmo_flux_sources:
                    flux_funcs[month][name][flux_source] = RectBivariateSpline(angles, e_grid, flux_data[flux_source])
                    # Note that e_grid is in GeV!

    with open("Fits/atmo_flux", "wb") as f:
        pickle.dump(flux_funcs, f)

    """
    Complete structure of flux_funcs should look like this:
    flux_funcs = {
        "January": {
            "_antinue": {
                "total": function,
                "pi": function,
                ...
            },
            "_nue": {
                ...
            }
        },
        "July": {
            ...
        }
    }
    """

def atmo_flux(E, theta, month, neutrino_name, flux_source):
    """Flux law for atmospheric neutrinos."""

    # MCEq defines theta=0 when neutrinos are coming down from the atmosphere, and theta>90 when they are up-going through the earth.
    # We have the opposite coordinate system, so we need to convert our angle first.
    theta = 180 * (1 - theta/np.pi)

    # Extracts the correct flux function associated with the last three parameters.
    flux_func = flux_funcs[month][neutrino_name][flux_source]

    # Computes the flux for a given angle and energy.
    flux = flux_func(theta, E)
    return flux[0][0]
