"""Collect a phase-model sweep into results.csv + results.png.
usage: collect_sector.py runs/phase_<design>_L1 [stroke_deg=3]"""
import sys, glob, json, csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

root = sys.argv[1].rstrip("/"); stroke = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
rows = sorted((json.load(open(f)) for f in glob.glob(f"{root}/*/summary.json")), key=lambda r: r["angle_deg"])
if not rows: sys.exit("no results")
a = np.array([r["angle_deg"] for r in rows]); T = np.array([r["torque_Nm"] for r in rows])

# best stroke-wide window (phase conducts for pitch/3 = 3 deg); interpolate torque on a fine grid
best = None
if a[-1] - a[0] >= stroke:
    fa = np.linspace(a[0], a[-1], 901); fT = np.interp(fa, a, T)
    for s0 in np.arange(a[0], a[-1] - stroke + 1e-9, 0.05):
        m = (fa >= s0 - 1e-9) & (fa <= s0 + stroke + 1e-9)
        avg = np.trapz(fT[m], fa[m]) / stroke
        if best is None or avg > best[1]: best = (s0, avg, fT[m].min())
psi = np.array([r["psi_phase_mWb"] for r in rows])
cols = ["angle_deg", "torque_Nm", "F_radial_N", "F_tangential_N", "psi_phase_mWb", "L_phase_mH",
        "phi_backiron_uWb", "B_backiron_T", "phi_tooth_uWb", "B_tooth_T", "B_core_max_T", "B_tooth_max_T", "ampere_turns_mesh", "n_tets"]
with open(f"{root}/results.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in rows: w.writerow([f"{r[c]:.5g}" for c in cols])

r0 = rows[0]
summary = dict(design=root.split("/")[-1], J_coil_A_mm2=r0["J_coil_A_mm2"], V_iron_cm3_per_pole=r0["V_core_cm3"] + r0["V_tooth_cm3"],
               T_peak_Nm=float(T.max()), angle_peak_deg=float(a[T.argmax()]),
               psi_aligned_over_unaligned=float(psi[0] / psi[-1]) if a[-1] >= 4.4 else None)
if best:
    summary.update(window_start_deg=float(best[0]), T_avg_window_Nm=float(best[1]), T_min_window_Nm=float(best[2]),
                   T_avg_per_cm3=float(best[1] / summary["V_iron_cm3_per_pole"]),
                   T_avg_phase_Nm=float(best[1] * 40))
json.dump(summary, open(f"{root}/design_summary.json", "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
ax[0].plot(a, T, "o-")
if best: ax[0].axvspan(best[0], best[0] + stroke, color="C2", alpha=.15, label=f"best {stroke:g}° window: avg {best[1]:.2f}, min {best[2]:.2f} N·m"); ax[0].legend()
ax[0].set(xlabel="rotor angle [deg]", ylabel="torque per pole [N·m]", title="Torque toward alignment (per pole; ×40 per phase)")
ax[1].plot(a, psi, "o-"); ax[1].set(xlabel="rotor angle [deg]", ylabel="ψ phase [mWb-turns]", title="Phase flux linkage (40 poles)")
for r, c in ((rows[0], "C0"), (rows[min(range(len(rows)), key=lambda i: abs(a[i] - 3))], "C1")):
    L = r["leg_sections"]
    ax[2].plot([x["r_mm"] for x in L], [x["B_T"] for x in L], "o-", color=c, label=f"{r['angle_deg']:g}° leg (mean B)")
    ax[2].axhline(r["B_backiron_T"], color=c, ls="--", lw=1, label=f"{r['angle_deg']:g}° back iron")
    ax[2].axhline(r["B_tooth_T"], color=c, ls=":", lw=1, label=f"{r['angle_deg']:g}° rotor tooth")
ax[2].set(xlabel="radius [mm] (back iron → pole tip)", ylabel="mean |B| over section [T]", title="Flux density along the stator leg"); ax[2].legend(fontsize=8)
for x in ax: x.grid(alpha=.3)
fig.suptitle(f"{summary['design']}  —  {r0['ampere_turns_target']:.0f} A-t, J = {r0['J_coil_A_mm2']:.1f} A/mm²")
fig.tight_layout(); fig.savefig(f"{root}/results.png", dpi=85)
print(open(f"{root}/results.csv").read()); print(json.dumps(summary, indent=1))
