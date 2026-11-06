"""Magnetic-circuit model for single-pole axial-flux SRM torque estimation.

Implements the nonlinear reluctance network described in the project notes.
"""

from __future__ import annotations

import numpy as np
from typing import Callable, Optional, Tuple


MU0 = 4 * np.pi * 1e-7  # H/m


class Material:
    """Simple material model. Supports constant mu or a B-H curve."""

    def __init__(
        self,
        mu_r: float = 1000.0,
        bh_curve: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    ):
        """
        Parameters
        ----------
        mu_r : float
            Constant relative permeability (used if bh_curve is None).
        bh_curve : tuple of arrays, optional
            (B, H) in SI units. If provided, overrides constant mu_r.
        """
        self.mu_r = mu_r
        self.bh_curve = bh_curve

    def mu_r_of_B(self, B: float) -> float:
        if self.bh_curve is None:
            return self.mu_r
        B_data, H_data = self.bh_curve
        B = np.clip(B, 0, B_data[-1])
        H = np.interp(B, B_data, H_data)
        if B < 1e-6:
            return self.mu_r
        return B / (MU0 * H)


class AxialPoleCircuit:
    """Single-pole axial-flux magnetic circuit with parallel rotor path.

    R is the average radius used for circumferential offset (d = R * theta).
    """

    def __init__(
        self,
        R: float,
        D: float,
        W: float,
        ell_s: float,
        ell_r: float,
        ell_a: float,
        N: int,
        material_s: Material,
        material_r: Material,
    ):
        self.R = R
        self.D = D
        self.W = W
        self.ell_s = ell_s
        self.ell_r = ell_r
        self.ell_a = ell_a
        self.N = N
        self.mat_s = material_s
        self.mat_r = material_r

    def _reluctances(self, phi: float, theta: float) -> Tuple[float, float, float]:
        """Return (R_s, R_a, R_r) at the given flux and position."""
        overlap = max(self.D - self.R * theta, 1e-9)
        non_overlap = max(self.R * theta, 1e-9)

        # Stator
        A_s = self.W * self.D
        B_s = abs(phi) / A_s if A_s > 0 else 0.0
        mu_s = self.mat_s.mu_r_of_B(B_s)
        R_s = self.ell_s / (mu_s * MU0 * A_s)

        # Air-gap (fixed)
        A_a = self.W * self.D
        R_a = self.ell_a / (MU0 * A_a)

        # Rotor parallel combination
        A_iron = self.W * overlap
        B_r = abs(phi) / A_iron if A_iron > 0 else 0.0
        mu_r = self.mat_r.mu_r_of_B(B_r)
        R_iron = self.ell_r / (mu_r * MU0 * A_iron)
        R_air = self.ell_r / (MU0 * self.W * non_overlap)
        R_r = 1.0 / (1.0 / R_iron + 1.0 / R_air)

        return R_s, R_a, R_r

    def solve_phi(self, i: float, theta: float, tol: float = 1e-8, max_iter: int = 60) -> float:
        """Solve for flux given current and position using bisection on residual."""
        if abs(i) < 1e-12:
            return 0.0

        target_mmf = self.N * i

        def residual(phi: float) -> float:
            R_s, R_a, R_r = self._reluctances(phi, theta)
            return (R_s + R_a + R_r) * phi - target_mmf

        # Bracket the root. Flux is positive for positive current.
        R_min = self.ell_a / (MU0 * self.W * self.D)
        phi_hi = target_mmf / R_min * 2.0
        phi_lo = 0.0

        for _ in range(10):
            if residual(phi_hi) > 0:
                break
            phi_hi *= 2.0

        for _ in range(max_iter):
            phi_mid = 0.5 * (phi_lo + phi_hi)
            r = residual(phi_mid)
            if abs(r) < tol * target_mmf or (phi_hi - phi_lo) < tol * phi_hi:
                return phi_mid
            if r > 0:
                phi_hi = phi_mid
            else:
                phi_lo = phi_mid
        return 0.5 * (phi_lo + phi_hi)

    def lambda_(self, i: float, theta: float) -> float:
        """Flux linkage."""
        return self.N * self.solve_phi(i, theta)

    def coenergy(self, i: float, theta: float, n_points: int = 40) -> float:
        """Numerical co-energy W_c = ∫ λ di' from 0 to i."""
        if abs(i) < 1e-12:
            return 0.0
        i_vals = np.linspace(0, i, n_points)
        lam_vals = np.array([self.lambda_(ii, theta) for ii in i_vals])
        return np.trapz(lam_vals, i_vals)

    def torque(self, i: float, theta: float, dtheta: float = 1e-4) -> float:
        """Torque via finite difference of co-energy at constant current."""
        Wc_plus = self.coenergy(i, theta + dtheta)
        Wc_minus = self.coenergy(i, theta - dtheta)
        return (Wc_plus - Wc_minus) / (2 * dtheta)


def linear_torque_analytic(
    i: float,
    theta: float,
    R: float,
    D: float,
    W: float,
    ell_s: float,
    ell_r: float,
    ell_a: float,
    N: int,
    mu_s: float,
    mu_r: float,
) -> float:
    """Closed-form torque for the constant-mu case (for validation)."""
    overlap = max(D - R * theta, 1e-9)
    non_overlap = max(R * theta, 1e-9)

    R_s = ell_s / (mu_s * MU0 * W * D)
    R_a = ell_a / (MU0 * W * D)
    R_iron = ell_r / (mu_r * MU0 * W * overlap)
    R_air = ell_r / (MU0 * W * non_overlap)
    R_r = 1.0 / (1.0 / R_iron + 1.0 / R_air)
    R_tot = R_s + R_a + R_r

    denom = mu_r * (D - R * theta) + R * theta
    dRr_dtheta = (ell_r / (MU0 * W)) * R * (mu_r - 1) / denom**2

    return -0.5 * (N * i / R_tot) ** 2 * dRr_dtheta
