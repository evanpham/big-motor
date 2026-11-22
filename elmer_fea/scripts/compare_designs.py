"""Compare phase-model design sweeps.
usage: compare_designs.py <out.png> runs/phase_A_L1 runs/phase_B_L1 ...
Runs collect_sector.py output (design_summary.json, summary.json per angle) must exist."""
import sys, glob, json, os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

out = sys.argv[1]; roots = [r.rstrip("/") for r in sys.argv[2:]]
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
hdr = f"{'design':26s} {'Tavg3°':>7s} {'Tmin3°':>7s} {'Tpeak':>6s} {'T/ph':>6s} {'ψa/ψu':>6s} {'Vfe':>5s} {'T/cm3':>6s} {'J':>5s} {'B0 tooth':>8s} {'B0 tip':>7s} {'B0 bi':>6s}"
print(hdr); print("-" * len(hdr))
for k, root in enumerate(roots):
    rows = sorted((json.load(open(f)) for f in glob.glob(f"{root}/*/summary.json")), key=lambda r: r["angle_deg"])
    S = json.load(open(f"{root}/design_summary.json"))
    lab = os.path.basename(root).replace("phase_", "").rsplit("_L", 1)[0]
    a = [r["angle_deg"] for r in rows]
    ax[0].plot(a, [r["torque_Nm"] for r in rows], "o-", ms=4, label=lab)
    ax[1].plot(a, [r["psi_phase_mWb"] for r in rows], "o-", ms=4, label=lab)
    L = rows[0]["leg_sections"]
    ax[2].plot([x["r_mm"] for x in L], [x["B_T"] for x in L], "-", color=f"C{k}", label=f"{lab} leg")
    ax[2].plot([L[-1]["r_mm"] + 0.8], [rows[0]["B_tooth_T"]], "s", color=f"C{k}", ms=8)
    ax[2].plot([L[0]["r_mm"] - 0.8], [rows[0]["B_backiron_T"]], "^", color=f"C{k}", ms=8)
    print(f"{lab:26s} {S.get('T_avg_window_Nm', np.nan):7.3f} {S.get('T_min_window_Nm', np.nan):7.3f} {S['T_peak_Nm']:6.3f} "
          f"{S.get('T_avg_phase_Nm', np.nan):6.1f} {S['psi_aligned_over_unaligned'] or np.nan:6.2f} {S['V_iron_cm3_per_pole']:5.2f} "
          f"{S.get('T_avg_per_cm3', np.nan):6.3f} {S['J_coil_A_mm2']:5.2f} {rows[0]['B_tooth_T']:8.2f} {L[-1]['B_T']:7.2f} {rows[0]['B_backiron_T']:6.2f}")
ax[0].set(xlabel="rotor angle [deg]", ylabel="torque per pole [N·m]", title="Torque toward alignment")
ax[1].set(xlabel="rotor angle [deg]", ylabel="ψ phase [mWb-turns]", title="Phase flux linkage")
ax[2].set(xlabel="radius [mm]   (▲ back iron · leg · ■ rotor tooth)", ylabel="mean |B| [T]", title="Aligned (0°) flux density along the path")
for x in ax: x.grid(alpha=.3); x.legend(fontsize=8)
fig.tight_layout(); fig.savefig(out, dpi=85)
print("\nTavg3°/Tmin3° = best 3° window average/minimum per pole [N·m]; T/ph = 40 × Tavg3°; Vfe = core + rotor tooth iron per pole [cm³];"
      "\nJ = gross coil current density [A/mm²]; B0 = mean |B| at 0° in the rotor tooth, the stator pole tip, and the back iron [T]")
