from NeutrinoFlux.imports import *

# Distance Functions
def x(theta):
    """Distance through the earth a particle must penetrate before reaching the detector."""
    return max(2 * R_E * np.cos(theta), 0)

def r(z, theta):
    """Returns the radius of the earth as a function of the angle theta and the penetration distance z = x'."""
    X = x(theta)
    r = (R_E**2 + (X - z)**2 - X * (X - z))**(0.5)
    return r

# Density Functions
avg_rho_earth = 5.51/(C.centi)**3 # mean mass density of the earth, in g/m^3

def rho_earth_average(r):
    """Returns the average mass density of the earth."""
    return avg_rho_earth

def rho_earth_full(r):
    """Returns the mass density at a point in the earth as a function of the radius. Taken from the Preliminary Earth Model."""
    y = r/R_E
    r /= C.kilo # convert to km
    
    rho = np.piecewise(
        r,
        [
            0 <= r < 1221.5,
            1221.5 <= r < 3480,
            3480 <= r < 5701,
            5701 <= r < 5771,
            5771 <= r < 5971,
            5971 <= r < 6151,
            6151 <= r < 6346.6,
            6346.6 <= r < 6356,
            6356 <= r < 6368,
            6368 <= r <= R_E/1e3
        ],
        [
            13.0885 - 8.8381*y*y,
            12.5815 - 1.2638*y - 3.6426*y*y - 5.5281*y*y*y,
            7.9565 - 6.4761*y + 5.5283*y*y - 3.0807*y*y*y,
            5.3197 - 1.4836*y,
            11.2494 - 8.0298*y,
            7.1089 - 3.8045*y,
            2.691 + 0.6924*y,
            2.9,
            2.6,
            1.02
        ]
    ) * 1/(C.centi)**3 # g/m**3
    
    return rho

def _rho_earth(z, theta):
    """Returns the nucleon number density as a function of penetration distance z and penetration angle theta."""
    return rho_earth_full(r(z, theta))

# Density Plotting
if __name__ == "__main__":
    for n in np.linspace(0, 0.9, 9):
        X = x(theta := np.pi * n)
        Z = np.linspace(0, X)
        rho = [_rho_earth(z, theta) for z in Z]
        plt.plot(Z, rho, label=f"$\\theta = {n:.1f}\pi$")

    plt.xlabel("Distance Penetrated (m)")
    plt.ylabel(r"Mass Density ($g/m^3$)")
    plt.title("Density Profile of the Earth")
    plt.legend()
    plt.show()

# Attenuation of the Earth
def attenuation_integral(theta):
    """Returns the integral of rho_earth described above."""
    return quad( _rho_earth, 0, x(theta), args=[theta] )[0]

angles = np.linspace(0, np.pi, 150)
attenuation = [attenuation_integral(angle) for angle in angles]
attenuation_function = InterpolatedUnivariateSpline(angles, attenuation, k=5) # Connects the dots

if __name__ == "__main__":
    with open("Fits/attenuation_function", "wb") as f:
        pickle.dump(attenuation_function, f)

    plt.plot(angles, [avg_rho_earth * x(angle) for angle in angles], label="Average", linestyle="dashed", color="green")
    plt.plot(angles, attenuation_function(angles), label="Actual Fit")
    plt.plot(angles, attenuation, linestyle="none", marker=".", markersize=3, label="Actual Points")

    plt.xlabel("$\\theta$ (rad)")
    plt.ylabel("$\\alpha(\\theta) := \int_0^{x(\\theta)}{\\rho(r(z, \\theta)) \\: d z}$ ($g/m^2$)")
    plt.title("Attenuation Parameter")
    plt.legend()
    plt.show()
