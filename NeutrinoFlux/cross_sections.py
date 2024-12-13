from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Callable, Self

import numpy as np

from . import constants as C

__all__ = ["cross_section_GR", "CrossSection"]

def cross_section_GR(E: C.Quantity) -> C.Quantity:
    """Gives the electron antineutrino cross section for GR events."""

    S_W = 2 * C.ureg.m_e * (E / C.ureg.c**2)/C.M_W**2
    sigma = C.sigma_0 * S_W/((1 - S_W)**2 + C.G_W**2)/(3 * C.R_W_mu)

    in_domain = (4 * C.ureg.PeV <= E) & (E <= 8 * C.ureg.PeV)
    return np.where(in_domain, sigma.to("cm**2"), 0)

@dataclass
class CrossSection:
    """Class containing the function and number of targets per mass for a given cross section."""

    name: str
    func: Callable[[C.Quantity], C.Quantity]
    target_concentration_earth: C.Quantity = C.nucleon_concentration_earth
    target_concentration_water: C.Quantity = C.nucleon_concentration_water

    cache = {}

    @classmethod
    def from_cache(cls, crsc_name: str) -> Self | None:
        return cls.cache.get(crsc_name, None)

    def __call__(self, E: C.Quantity | float, units: bool = True) -> C.Quantity | float:
        """Evaluates the cross section for a given energy."""

        if not units:
            return self.func(E * C.ureg.GeV).m
        return self.func(E)

    def times_conc(self, E: C.Quantity, target: str) -> C.Quantity:
        crsc = self(E)
        if target == "earth":
            return crsc * self.target_concentration_earth
        if target in ("ice", "water"):
            return crsc * self.target_concentration_water
        raise NotImplementedError

    def plot(self, ax, E: C.Quantity | None = None):
        """Plots sigma(E) vs E."""

        if E is None:
            E = np.logspace(1, 12, 111) * C.ureg.GeV
        sigma = self(E)

        ax.set_xscale("log")
        ax.set_yscale("log")

        plot, = ax.plot(E, sigma, label=self.name)
        return plot

@C.ureg.wraps("cm**2", ["GeV", "cm**2"])
def interpolate_sigma(E: float, sigma_p: np.ndarray) -> float:
    """Interpolates default ISO cross section data."""

    log_E_p = np.linspace(1, 12, 111)
    log_E = np.log10(E)
    log_sigma_p = np.log10(sigma_p)

    log_sigma = np.interp(log_E, log_E_p, log_sigma_p, 0, 0)
    return 10 ** log_sigma

for file in (Path(__file__).parent/"data/cross_section").iterdir():
    name = file.stem.removeprefix("total_").removesuffix("_iso_NLO_HERAPDF1.5NLO_EIG").lower()
    sigma_ = 1e-27 * C.ureg.cm**2 * np.loadtxt(file, skiprows=1, unpack=True)
    func = partial(interpolate_sigma, sigma_p=sigma_)

    CrossSection.cache[name] = CrossSection(name, func)

CrossSection.cache["gr"] = CrossSection("gr", cross_section_GR,
                                        C.electron_concentration_earth,
                                        C.electron_concentration_water)
