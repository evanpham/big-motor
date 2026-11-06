"""Example: torque vs rotor position for a few currents."""

import numpy as np
import matplotlib.pyplot as plt
from mag_analysis import AxialPoleCircuit, Material

# Geometry (SI units)
geom = dict(
    R=0.15,       # outer radius [m]
    D=0.04,       # pole width [m]
    W=0.03,       # axial length [m]
    ell_s=0.08,
    ell_r=0.06,
    ell_a=0.001,  # air-gap
    N=100,
)

mat = Material(mu_r=2000)  # linear for this demo

circuit = AxialPoleCircuit(**geom, material_s=mat, material_r=mat)

thetas = np.linspace(0.001, 0.25, 40)  # rad
currents = [5, 10, 20]

plt.figure(figsize=(8, 5))
for i in currents:
    torques = [circuit.torque(i, th) for th in thetas]
    plt.plot(np.degrees(thetas), torques, label=f"{i} A")

plt.xlabel("Rotor position [deg]")
plt.ylabel("Torque [N·m]")
plt.title("Single-pole axial SRM – linear model")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("torque_vs_theta.png", dpi=150)
print("Saved torque_vs_theta.png")
