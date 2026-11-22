"""Post-process one run of the multi-pole phase model (sector.geo / case_sector.sif).
Frame: motor axis = z, middle stator pole on +x, rotor tooth rotated by +theta toward +y.
usage: post_sector.py <run dir> <theta deg> <turns> <amps> <coil area m2>"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readvtu, probe

run = sys.argv[1]; th = np.radians(float(sys.argv[2])); N = float(sys.argv[3]); I = float(sys.argv[4]); Acoil = float(sys.argv[5])
v = readvtu.read_results(f"{run}/results")
X = v["points"]; P = v["point"]; T = v["types"]
tet = np.array([c for c, t in zip(v["cells"], T) if t == 10]); g = v["cell"]["GeometryIds"][T == 10]
B = P["magnetic flux density e"]; J = P["current density e"]
a, b, c, d = (X[tet[:, k]] for k in range(4)); Vt = np.abs(np.einsum('ij,ij->i', b-a, np.cross(c-a, d-a))) / 6
pr = {k: probe.Probe(X, tet[g == k]) for k in (1, 2, 3)}

def surf(body, pts, nrm, dA):
    """flux of the field through sample points (one area element each) inside a body; returns flux, area"""
    f, ok = pr[body].interp(B if body != 2 else J, pts)
    return float(np.nansum(np.einsum('ij,ij->i', np.nan_to_num(f), nrm))) * dA, float(ok.sum()) * dA

n = 220
# coil current through the half-plane y=0 (azimuthal direction there is +y)
u = np.linspace(80e-3, 100e-3, n); w = np.linspace(-15e-3, 15e-3, n); U, W = np.meshgrid(u, w)
Icut, _ = surf(2, np.c_[U.ravel(), np.zeros(U.size), W.ravel()], np.tile([0, 1, 0], (U.size, 1)), (u[1]-u[0])*(w[1]-w[0]))

# back iron of the middle core: axial flux through z=0 (only core cells count)
u = np.linspace(70e-3, 100e-3, n); w = np.linspace(-7e-3, 7e-3, n); U, W = np.meshgrid(u, w)
phi_bi, A_bi = surf(1, np.c_[U.ravel(), W.ravel(), np.zeros(U.size)], np.tile([0, 0, 1], (U.size, 1)), (u[1]-u[0])*(w[1]-w[0]))
phi_bi = abs(phi_bi)
_, okb = pr[1].interp(B, np.c_[U.ravel(), W.ravel(), np.zeros(U.size)])
r_bi_out = U.ravel()[okb].max()          # outer radius of the back iron

# one leg (z > 0) of the middle core: radial flux through cylinders r = R
legs = []
for R in np.r_[np.arange(np.ceil((r_bi_out + 0.3e-3)*2e3)/2e3, 99.5e-3, 0.5e-3), 99.6e-3]:
    ph = np.linspace(-5.5, 5.5, n) * np.pi/180; zz = np.linspace(0, 30e-3, n); PH, ZZ = np.meshgrid(ph, zz)
    pts = np.c_[R*np.cos(PH.ravel()), R*np.sin(PH.ravel()), ZZ.ravel()]
    nr = np.c_[np.cos(PH.ravel()), np.sin(PH.ravel()), np.zeros(PH.size)]
    f, A = surf(1, pts, nr, R*(ph[1]-ph[0])*(zz[1]-zz[0]))
    if A > 0: legs.append(dict(r_mm=round(R*1e3, 2), area_mm2=A*1e6, flux_uWb=abs(f)*1e6, B_T=abs(f)/A))

# middle rotor tooth: axial flux through z=0, sampled in the tooth's own frame
ct, st = np.cos(th), np.sin(th)
# radial range must cover the deepest rotor tooth (was 100-112 mm, which cut off teeth deeper than 12 mm)
u = np.linspace(100e-3, 130e-3, 2*n); w = np.linspace(-4e-3, 4e-3, n); U, W = np.meshgrid(u, w)
pts = np.c_[U.ravel()*ct - W.ravel()*st, U.ravel()*st + W.ravel()*ct, np.zeros(U.size)]
phi_rt, A_rt = surf(3, pts, np.tile([0, 0, 1], (U.size, 1)), (u[1]-u[0])*(w[1]-w[0]))
phi_rt = abs(phi_rt)

sf = f"{run}/results/scalars.dat"
s = np.loadtxt(sf if os.path.exists(sf) else sf + ".0", ndmin=1)
# columns: eddy power, energy, F stator xyz, T stator, F rotor xyz, T rotor
Tst, Trt = s[5], s[9]
Fx, Fy = s[6], s[7]
Bm = np.linalg.norm(B[tet].mean(1), axis=1)
leg_min = min(legs, key=lambda L: L["area_mm2"]) if legs else {}
res = dict(
    angle_deg=float(sys.argv[2]), n_tets=int(len(tet)),
    ampere_turns_target=N*I, ampere_turns_mesh=float(abs(Icut)), coil_area_mm2=Acoil*1e6, J_coil_A_mm2=N*I/Acoil/1e6,
    torque_Nm=float(0.5*(Tst - Trt)),            # per pole, positive = pulls the rotor back toward 0 deg
    torque_stator_Nm=float(Tst), torque_rotor_Nm=float(Trt),
    F_radial_N=float(Fx*ct + Fy*st), F_tangential_N=float(-Fx*st + Fy*ct),   # on the middle rotor tooth
    r_backiron_out_mm=float(r_bi_out*1e3), phi_backiron_uWb=phi_bi*1e6, A_backiron_mm2=A_bi*1e6, B_backiron_T=phi_bi/A_bi,
    psi_pole_mWb=N*phi_bi*1e3, psi_phase_mWb=40*N*phi_bi*1e3, L_phase_mH=40*N*phi_bi/I*1e3,
    phi_tooth_uWb=phi_rt*1e6, A_tooth_mm2=A_rt*1e6, B_tooth_T=phi_rt/A_rt,
    B_core_max_T=float(Bm[g == 1].max()), B_tooth_max_T=float(Bm[g == 3].max()),
    V_core_cm3=float(Vt[g == 1].sum()*1e6), V_tooth_cm3=float(Vt[g == 3].sum()*1e6),
    leg_sections=legs)
json.dump(res, open(f"{run}/summary.json", "w"), indent=1)
for k, val in res.items():
    if k != "leg_sections": print(f"  {k:22s} {val:.5g}")
print("  leg (z>0), radial flux through r = const:")
for L in legs: print(f"    r={L['r_mm']:6.2f} mm  A={L['area_mm2']:6.1f} mm2  flux={L['flux_uWb']:6.1f} uWb  B={L['B_T']:.2f} T")
