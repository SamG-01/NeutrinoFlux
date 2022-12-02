# NeutrinoFlux
Computes the yearly rate of neutrino events in the IceCube Neutrino Observatory at the south pole.

## Basic Usage
As a quick example, the code below computes the yearly rate of electron neutrinos from astrophysical sources, from all zenith angles (relative to the north pole at 0 radians, with the south pole at the origin), for energies between 10<sup>13</sup> and 10<sup>21</sup> eV.

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

If I wanted to compute only the rate of Glashow Resonance events for electron antineutrinos within the same range, but with a spectral index of 2 (instead of the default 2.53) and a normalization of 1 (instead of the default 1.66), I would instead write:
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
print(rate)
```

If I instead consider atmospheric tau neutrino flux from pion decay in July with energies up to 10<sup>16</sup> eV, I would write:
```python
rate = event_rate(
    nu=nu["nu_tau"],
    flux_type="atmo",
    E_bounds=(1e13, 1e16),
    theta_bounds=full_theta_bounds,
    diff_flux_kwargs={
        "month": "July",
        "atmo_flux_source": "pr"
    }
)
print(rate)
```

## Advanced Usage
TODO

## Limitations
TODO
