from NeutrinoFlux.imports import *

# CrossSection class
class CrossSection():
    """Class containing the functions, bounds, and number of targets per mass for a given cross section."""
    def __init__(self, name, anti, func, E_domain, targets_per_gram_earth, targets_per_gram_water) -> None:
        # basic properties: name, and whether it applies to neutrinos or antineutrinos
        self.name = name
        self.anti = anti
        
        # sigma(E) function and its domain (domain is all energies if True). sigma(E) = 0 outside of this domain
        self.func = func # input is E in eV, output is sigma in m^2
        self.E_domain = E_domain

        # targets per gram
        self.targets_per_gram_earth = targets_per_gram_earth
        self.targets_per_gram_water = targets_per_gram_water

    def eval(self, E):
        """Evaluates the cross section at given bounds."""
        if self.E_domain == True:
            return self.func(E)

        E_a, E_b = self.E_domain
        if E_a <= E <= E_b:
            return self.func(E)
        else: return 0

# GR Contribution
GR_bounds = GR_a, GR_b = (4e15, 8e15)

def cross_section_GR(E):
    """Gives the electron antineutrino cross section for GR events."""
    GF2 = 1.3604656e-10 * (1e9)**(-4) # eV**(-4)
    MW2 = 6460.783641 * (1e9)**2 # eV**2
    GeV2_MBARN = 0.3893796623 * (1e9)**2 * 1e-27 * (1e-2)**2 # eV**2 * m**2
    GW2 = 6.935717E-4 # unitless
    RWmu =  0.1057 # unitless

    crs0 = GF2*MW2/np.pi*GeV2_MBARN
    m_electron = 5.10999E-4 * 1e9 # eV
    SW = 2*m_electron*E/MW2
    sigma = crs0*SW/( (1 - SW)*(1 - SW) + GW2)/RWmu/3
    return sigma

# GR Plotting
if __name__ == "__main__":
    E_res = 6.3e15 #* eV
    GR_range = np.linspace(4e15, 8e15)
    plt.plot(GR_range, cross_section_GR(GR_range), label="GR")
    plt.axvline(x=E_res, linestyle="dashed", label="Resonance Energy")

    #plt.xscale('log')
    plt.yscale('log')

    plt.title("GR Cross Section")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Cross Section (m$^2$)")

    plt.legend()
    plt.show()

    print("Resonant Cross Section:", max(cross_section_GR(GR_range)))

# Performs fitting for the NC and CC cross sections.
files = [path + "/Cross Section Data/total_nu" + name + "_iso_NLO_HERAPDF1.5NLO_EIG.dat" for name in ("_CC", "_NC", "bar_NC", "bar_CC")]

nu_cc, nu_nc, nubar_nc, nubar_cc = 1e-27 * (C.centi)**2 * np.array([
    np.loadtxt(file, skiprows=1, unpack=True) for file in files
])

E_range = np.logspace(1, 12, 111) * C.giga

cross_section_data = {
    False: [nu_cc, nu_nc, "Neutrino"],
    True: [nubar_cc, nubar_nc, "Antineutrino"]
} # argument: anti

cross_sections = {}

for anti in [False, True]:
    cc, nc, name = cross_section_data[anti]

    cc_fit = InterpolatedUnivariateSpline(E_range, cc, k=5)
    nc_fit = InterpolatedUnivariateSpline(E_range, nc, k=5)
    tot_fit = InterpolatedUnivariateSpline(E_range, cc + nc, k=5)

    cross_sections[anti] = {
        "nc": CrossSection("nc", anti, nc_fit, True, nucleons_per_gram_earth, nucleons_per_gram_water),
        "cc": CrossSection("cc", anti, cc_fit, True, nucleons_per_gram_earth, nucleons_per_gram_water)
    }

    cross_sections["GR"] = CrossSection("GR", anti, cross_section_GR, GR_bounds, electrons_per_gram_earth, electrons_per_gram_water)

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
