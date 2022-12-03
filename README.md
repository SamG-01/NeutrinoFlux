# NeutrinoFlux
Predicts the yearly rate of neutrino events in the IceCube Neutrino Observatory at the south pole. For this purpose, this package also includes modules with functions for neutrino cross sections, differential flux from atmospheric and astrophysical sources, and attenuation through the earth.

To allow for compatibility with other packages, this package uses classes for `Neutrino`s and `CrossSection`s, whose default attributes and methods can be replaced with other functions as desired.

## Basic Usage
As a quick example, the code below computes the yearly rate of electron neutrinos from astrophysical sources, from all zenith angles, for energies between $10^{13}$ and $10^{21}$ eV.

```python
import numpy as np
from NeutrinoFlux.event_rate import event_rate
from NeutrinoFlux.neutrinos import default_neutrinos as nu

full_E_bounds = (1e13, 1e21)
full_theta_bounds = (0, np.pi)

rate = event_rate(
    nu=nu["nu_e"],
    flux_type="astro",
    E_bounds=full_E_bounds,
    theta_bounds=full_theta_bounds
)
print(rate)
```

If I wanted to compute only the rate of Glashow Resonance events for electron antineutrinos within the same range, but with a spectral index of $\gamma = 2$ (instead of the default $2.53$) and a normalization of $\phi_{\text{astro}} = 1$ (instead of the default $1.66$), I would instead write:
```python
rate = event_rate(
    nu=nu["nubar_e"],
    flux_type="astro",
    E_bounds=full_E_bounds,
    theta_bounds=full_theta_bounds,
    diff_flux_kwargs={
        "gamma": 2,
        "phi_astro": 1,
    },
    GR_only=True
)
```

If I instead consider atmospheric tau neutrino flux from pion decay in July with energies only up to $10^{16}$ eV, I would write:
```python
rate = event_rate(
    nu=nu["nu_tau"],
    flux_type="atmo",
    E_bounds=(1e13, 1e16),
    theta_bounds=full_theta_bounds,
    diff_flux_kwargs={
        "month": "July",
        "atmo_flux_source": "pi"
    }
)
```

## Background
TODO

- link paper?

## Advanced Usage
TODO

- zenith angle: relative to the north pole at $\theta = 0$ radians, with the south pole at the origin)
- discuss CrossSection, Neutrino classes, making your own modifications

## Limitations
TODO

- fits have a minimum of $10^{13}$ eV: that's just my default models tho
- tau and mu flux may be inaccurate due to other effects

## Compatibility
While this package provides a set of models by default, other models can be incorporated by modifying a `Neutrino` object's attributes.

For example, the [`nuFATE`](https://github.com/aaronvincent/nuFATE) can be used to compute the attenuation of initial flux for a given flavor, including $\nu_{\tau}$ regeneration effects. To use this alongside `NeutrinoFlux`, all it would take is taking a `Neutrino`-class object and replacing its default `attenuation` method with one that uses the `nuFATE` attenuation model.
