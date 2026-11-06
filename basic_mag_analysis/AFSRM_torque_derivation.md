# Axial-Flux SRM Single-Pole Torque Derivation
## (Magnetic-circuit / co-energy approach)

### Geometry and notation
- \(R\): **average** radius of the pole arc, \(R = (R_\text{outer} + R_\text{inner})/2\)  
  (used so that circumferential offset \(d = R\theta\); this is more accurate than using the outer radius alone)
- \(D\): pole width (aligned overlap dimension)
- \(W\): axial length (or stack width) of the pole
- \(\ell_s, \ell_r, \ell_a\): magnetic path lengths in stator iron, rotor iron, and main air-gap
- \(\mu_s, \mu_r\): relative permeabilities of stator and rotor iron
- \(\mu_0\): permeability of free space
- \(N\): number of turns
- \(\theta\): rotor position (0 = aligned)
- \(\mathcal{R}\): total reluctance of the magnetic circuit (double-struck R in the original notes)

Percent alignment is defined as
\[
\%\text{alignment} = 100 \times \left(1 - \frac{R\theta}{D}\right).
\]

### Assumptions
1. All flux is guided by the iron; leakage and fringing outside the simple circuit are neglected.
2. Flux is parallel to the core path (no radial or axial spreading inside the iron).
3. The rotor section is modelled as the **parallel** combination of:
   - iron path over the remaining overlap width \((D - R\theta)\)
   - air path over the misaligned width \(R\theta\)
4. Main stator–rotor air-gap reluctance is treated as constant (fixed gap length and nominal area).
5. Iron permeabilities are in general nonlinear: \(\mu = \mu(B)\) taken from a material \(B\)–\(H\) (or magnetization) curve. Consequently the iron reluctances (and the total circuit reluctance) depend on the operating point.
6. Single pole pair; mutual coupling to other poles is ignored.

### Reluctance network

Iron reluctances are operating-point dependent. Given a material curve, one typically works with \(B\) (or \(H\)) in each iron segment:

- Guess or iterate on the flux \(\phi\) (or \(\lambda\)) until the MMF drop around the loop equals \(Ni\).
- Local \(B = \phi / A\) in each iron segment determines \(\mu(B)\) or \(H(B)\) from the datasheet.
- Segment reluctance \(\mathcal{R} = \ell H / \phi\) (or \(\ell /(\mu(B)\mu_0 A)\)).

**Stator iron** (nonlinear, depends on local \(B_s\)):
\[
\mathcal{R}_s(B_s) = \frac{\ell_s}{\mu_s(B_s)\,\mu_0 W D}
\]

**Main air-gap** (linear, fixed geometry):
\[
\mathcal{R}_a = \frac{\ell_a}{\mu_0 W D}
\]

**Rotor section** (parallel iron + air; iron part nonlinear and position-dependent):
\[
\mathcal{R}_{\text{iron}}(B_r,\theta) = \frac{\ell_r}{\mu_r(B_r)\,\mu_0 W(D - R\theta)}, \qquad
\mathcal{R}_{\text{air}}(\theta) = \frac{\ell_r}{\mu_0 W\cdot R\theta}
\]
\[
\frac{1}{\mathcal{R}_r} = \frac{1}{\mathcal{R}_{\text{iron}}} + \frac{1}{\mathcal{R}_{\text{air}}}
\quad\Rightarrow\quad
\mathcal{R}_r(\phi,\theta) = \frac{\ell_r}{\mu_0 W\bigl[\mu_r(B_r)(D - R\theta) + R\theta\bigr]}
\]
(with \(B_r = \phi / [W(D-R\theta)]\) on the iron branch).

**Total reluctance** (implicit function of current via the operating-point flux):
\[
\mathcal{R}(i,\theta)\quad\text{or equivalently}\quad\mathcal{R}(\lambda,\theta)
\]
obtained by solving the nonlinear magnetic circuit so that
\[
Ni = \mathcal{R}_s\phi + \mathcal{R}_a\phi + \mathcal{R}_r\phi.
\]

### Flux linkage
By definition
\[
\lambda(i,\theta) = N\phi(i,\theta),
\]
where \(\phi\) is the solution of the nonlinear circuit above.  
Equivalently, once the operating-point reluctance is known,
\[
\lambda = \frac{N^2 i}{\mathcal{R}(i,\theta)}.
\]
(The second form is convenient after the circuit has been solved; it is not an independent closed-form expression when \(\mathcal{R}\) depends on \(i\).)

### Co-energy (general)
Co-energy is defined by the integral at constant rotor position:
\[
W_c(i,\theta) = \int_0^i \lambda(i',\theta)\,di'.
\]
This expression is valid for both linear and nonlinear media.  
(The familiar \(\frac12\lambda i\) or \(\frac{(Ni)^2}{2\mathcal{R}}\) formulas hold **only** when \(\mathcal{R}\) is independent of current.)

### Torque (virtual work, constant current)
\[
T_e = \left.\frac{\partial W_c}{\partial\theta}\right|_{i=\text{const}}.
\]
In practice this is evaluated numerically:
1. For a grid of currents and positions, solve the nonlinear circuit → obtain \(\lambda(i,\theta)\).
2. Integrate with respect to current (e.g. trapezoidal rule) to obtain \(W_c(i,\theta)\).
3. Differentiate \(W_c\) with respect to \(\theta\) at fixed \(i\) (finite differences or spline derivative).

### Linear special case (constant \(\mu\))
When permeabilities are constant the circuit is linear, \(\mathcal{R}=\mathcal{R}(\theta)\) only, and the expressions collapse to
\[
\lambda = \frac{N^2 i}{\mathcal{R}(\theta)}, \qquad
W_c = \frac{(Ni)^2}{2\mathcal{R}(\theta)},
\]
\[
T_e = -\frac12\left(\frac{Ni}{\mathcal{R}}\right)^2\frac{\partial\mathcal{R}}{\partial\theta}.
\]
The analytic rearranged form given in the original notes applies in this limit.

### Implementation notes for code
- Guard the domain: \(0 \le R\theta < D\) (or switch to a different expression past full misalignment).
- Nonlinear solution: for each \((i,\theta)\) iterate on \(\phi\) (or \(\lambda\)) until MMF balance is satisfied, using the material \(B\)–\(H\) curve for the iron segments. A robust root-finding method (bisection or Newton) on the MMF residual is preferred over simple fixed-point iteration once the material saturates.
- Prefer working from the \(B\)–\(H\) datasheet directly rather than a pre-fitted \(\mu(B)\); convert to \(\mu = B/(\mu_0 H)\) only if needed.
- Co-energy must be obtained by numerical integration of \(\lambda(i)\) whenever the material is nonlinear.
- At exact alignment (\(\theta = 0\)) torque is identically zero by symmetry; numerical finite differences should special-case this point.
- Expected limitations (for later FEA comparison):
  - Neglect of 3-D leakage/fringing underestimates \(\lambda\) in the unaligned position → overestimates torque at large \(\theta\).
  - The average-radius approximation is better than the outer-radius approximation; an integrated-arc treatment is still more accurate.
  - The simple parallel iron/air split is only a first-order representation of the position-dependent flux paths.

### Quick validation checks
- At \(\theta=0\) (fully aligned) the air branch disappears and \(\mathcal{R}_r\to\ell_r/(\mu_r(B_r)\mu_0 WD)\).
- As \(R\theta\to D\) the iron branch reluctance \(\to\infty\) and \(\mathcal{R}_r\to\ell_r/(\mu_0 W R\theta)\).
- In the linear limit the analytic torque formula is recovered.
- \(T_e\) must be an odd function about the aligned position and zero at both fully aligned and fully unaligned extremes for a symmetric pole.
