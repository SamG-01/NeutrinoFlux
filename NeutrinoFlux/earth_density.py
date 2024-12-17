import numpy as np

from . import constants as C

# Distance Functions


def penetration_length(theta: C.Quantity) -> C.Quantity:
    """Returns the distance through the earth a particle must
    penetrate before reaching the detector."""

    return np.clip(2 * C.R_E * np.cos(theta), 0, None)


def r_E(z: C.Quantity, theta: C.Quantity,
        x: C.Quantity | None = None) -> C.Quantity:
    """Returns distance from the center of the earth as a function
    of the angle theta and the penetration distance z."""

    if x is None:
        x = penetration_length(theta)
    return ((C.R_E**2 + (x - z)**2 - x * (x - z))**0.5).to(z.u)


@C.ureg.wraps("g/cm**3", "km")
def rho_earth(r: float) -> float:
    """Returns the mass density at a point in the earth as a
    function of the radius. Taken from the Preliminary Earth Model."""

    y = r/C.R_E.m

    d = np.select(
        [
            (0 <= r) & (r < 1221.5),
            (1221.5 <= r) & (r < 3480),
            (3480 <= r) & (r < 5701),
            (5701 <= r) & (r < 5771),
            (5771 <= r) & (r < 5971),
            (5971 <= r) & (r < 6151),
            (6151 <= r) & (r < 6346.6),
            (6346.6 <= r) & (r < 6356),
            (6356 <= r) & (r < 6368),
            (6368 <= r) & (r <= C.R_E.m)
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
        ],
        0
    )

    return d


def attenuation_parameter(theta: C.Quantity, N: int = 1000) -> C.Quantity:
    x = penetration_length(theta)
    z = np.linspace(0, x, N)
    r = r_E(z, theta, x)
    return np.trapz(rho_earth(r), z, axis=0)


def attenuation_plot(ax) -> tuple:
    """Plots alpha(theta) vs theta."""

    theta = np.linspace(0, np.pi, 150) * C.ureg.rad
    average, = ax.plot(theta, C.rho_earth_avg * penetration_length(theta),
                       label="Average", linestyle="dashed")
    actual, = ax.plot(theta, attenuation_parameter(theta), label="Actual")

    return average, actual
