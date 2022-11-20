from NeutrinoFlux.imports import *

# GR Contribution
GR_bounds = GR_a, GR_b = (4*C.peta, 8*C.peta)

def _sigma_GR(E):
    """Gives the electron antineutrino cross section for GR events. Source: PDG 2021."""
    #E /= C.giga # Convert E from eV to GeV
    GF2 = 1.3604656e-10 * (C.giga)**(-4) # eV**(-4)
    #    MW2 = 6467.858929
    #PDG 2021
    MW2 = 6460.783641 * (C.giga)**2 # eV**2
    GeV2_MBARN = 0.3893796623 * (C.giga)**2 * 1e-27 * (C.centi)**2 # eV**2 * m**2
    GW2 = 6.935717E-4 # unitless
    RWmu =  0.1057 # unitless

    crs0 = GF2*MW2/np.pi*GeV2_MBARN
    m_electron = 5.10999E-4 * C.giga # eV
    SW = 2*m_electron*E/MW2
    sigma = crs0*SW/( (1 - SW)*(1 - SW) + GW2)/RWmu/3
    return sigma

def cross_section_GR(E):
    return _sigma_GR(E) if GR_a <= E <= GR_b else 0

# GR Plotting
if __name__ == "__main__":
    E_res = 6.3e15 #* eV
    GR_range = np.linspace(0, 10e15)
    plt.plot(GR_range, _sigma_GR(GR_range), label="GR")
    plt.axvline(x=E_res, linestyle="dashed", label="Resonance Energy")

    #plt.xscale('log')
    plt.yscale('log')

    plt.title("GR Cross Section")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Cross Section (m$^2$)")

    plt.legend()
    plt.show()

    print("Resonant Cross Section:", max(_sigma_GR(GR_range)))

# NC and CC Contribution
files = [path + "/Cross Section Data/total_nu" + name + "_iso_NLO_HERAPDF1.5NLO_EIG.dat" for name in ("_CC", "_NC", "bar_NC", "bar_CC")]

nu_cc, nu_nc, nubar_nc, nubar_cc = 1e-27 * (C.centi)**2 * np.array([
    np.loadtxt(file, skiprows=1, unpack=True) for file in files
])

E_range = np.logspace(1, 12, 111) * C.giga

cross_section_data = {
    False: [nu_cc, nu_nc, "Neutrino"],
    True: [nubar_cc, nubar_nc, "Antineutrino"]
} # argument: anti

cross_sections = {
    
}

for anti in [False, True]:
    cc, nc, name = cross_section_data[anti]

    cc_fit = InterpolatedUnivariateSpline(E_range, cc, k=5)
    nc_fit = InterpolatedUnivariateSpline(E_range, nc, k=5)
    tot_fit = InterpolatedUnivariateSpline(E_range, cc + nc, k=5)

    cross_sections[anti] = {
        "nc": nc_fit,
        "cc": cc_fit,
        "tot": tot_fit
    }

    # Cross Section Plotting
    if __name__ == "__main__":
        plt.plot(E_range, cc_fit(E_range), label="CC")
        plt.plot(E_range, nc_fit(E_range), label="NC")
        plt.plot(E_range, tot_fit(E_range), label="TOT")

        plt.plot(E_range, cc, label="CC Points", marker=".", ls="none", markersize=1)
        plt.plot(E_range, nc, label="NC Points", marker=".", ls="none", markersize=1)
        plt.plot(E_range, cc + nc, label="TOT Points", marker=".", ls="none", markersize=1)

        plt.xscale('log')
        plt.yscale('log')

        plt.title(name + " Cross Sections")
        plt.xlabel("Energy (eV)")
        plt.ylabel("Cross Section (m$^2$)")
        
        plt.legend()
        plt.show()

if __name__ == "__main__":
    with open("Fits/cross_sections", "wb") as f:
        pickle.dump(cross_sections, f)
