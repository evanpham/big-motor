# mag_analysis

Magnetic-circuit library for estimating torque of a single-pole axial-flux switched reluctance machine.

Based on the reluctance network with a parallel iron/air path in the rotor section. Supports both constant-permeability and nonlinear $B$–$H$ materials.

## Install

```bash
pip install -e .
```

## Quick start

```python
from mag_analysis import AxialPoleCircuit, Material

mat = Material(mu_r=1500)
circuit = AxialPoleCircuit(
    R=0.15, D=0.04, W=0.03,
    ell_s=0.08, ell_r=0.06, ell_a=0.001,
    N=80,
    material_s=mat,
    material_r=mat,
)

lam = circuit.lambda_(10.0, theta=0.05)
T = circuit.torque(10.0, theta=0.05)
print(lam, T)
```

See `examples/` for plotting scripts.
