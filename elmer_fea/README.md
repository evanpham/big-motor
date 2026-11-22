# elmer_fea: 3D FEA of one SRM phase in Elmer

3D magnetostatic finite-element analysis of one phase of the 40-pole ring-coil machine, used to compare design variants: rotor tooth depth, back-iron thickness and phase current. Each run computes static torque, flux linkage and flux density at one rotor angle. A sweep over 0–4.5° gives the torque curve for one design.

The machine, as modelled: one ring coil (45 turns) threads 40 stator C-cores. The outer rotor has I-bar teeth on the same 9° pitch. The gap is at 200 mm diameter and 0.5 mm wide. The baseline current is 36 A.

## Model

- **3-pole segment.** The model contains 3 stator C-cores (9° apart), the ring coil cut to the segment, every rotor tooth inside it, and air out to r = 40–160 mm and z = ±80 mm. Torque is taken on the **middle rotor tooth only** (Maxwell stress), so the neighbouring poles supply the correct surroundings, including their parasitic pull. The cut faces are flux-parallel (n × A = 0).
- **Validation.** The middle-pole torque from the 3-pole model is within 0.4% of a 5-pole model. A true periodic 1-pole sector did not converge with edge elements in this Elmer build (9.0).
- **Mesh.** Gmsh (HXT), with sizes set per region. `L1` (all results here): gap 0.5 mm, iron 1.2 mm, coil 2.5 mm. L1 is within about 1% of the finer `L2` mesh, so differences below about 1% between designs are not resolved.
- **Solver.** Elmer `WhitneyAVSolver` (edge elements), nonlinear M-19 B-H curve with Newton iterations, 4 MPI processes. About 4 minutes per angle on a 4-core i7.
- **Coil.** The CAD coil grazes the core fillets. The geometry script trims the overlap (about 40 mm³) and reports it. The current density is computed from the trimmed coil area (214.85 mm² → 7.54 A/mm² at 36 A).

## Running

Requires `gmsh`, Elmer (`ElmerSolver_mpi`, `ElmerGrid`), `mpirun` and Python 3 with numpy, scipy and matplotlib.

```bash
cd elmer_fea
scripts/sweep_sector.sh designs/13_10_5.6_10_3.5_2_3.zip            # 0, 0.5 ... 4.5 deg -> runs/phase_<design>_L1/
scripts/sweep_sector.sh designs/13_10_5.6_10_3.5_2_3.zip 0.5 2.5     # chosen angles only
AMPS=50 scripts/run_sector.sh designs/13_10_5.6_10_3.5_2_3 2.5 phase_13_..._A50/02.5deg   # one run, other current
python3 scripts/collect_sector.py runs/phase_<design>_L1             # results.csv + design_summary.json
python3 scripts/compare_designs.py out.png runs/phase_A_L1 runs/phase_B_L1 ...
python3 -I scripts/leakage_lines.py runs/phase_<design>_L1/00.0deg out.png   # leakage field-line figure
```

Environment options: `LEVEL=L1|L2`, `NP` (MPI processes), `AMPS`, `TURNS`, and `MIN_FRAG` (see the caveats below). Sweeps skip angles that already have a `summary.json`, so an interrupted sweep resumes where it stopped.

**Design names** follow `designs/namingConvention.txt`: `<rotor_tooth_depth>_<back_iron_depth>_<stator_tooth_width>_<stator_leg_height>_<BackIronAdditionalThk>_<taper_step_length>_<inner_fillet_rad>` (mm). Each zip holds the Onshape STEP export of the stator C-core and the rotor I-bar, at 0° (aligned). The ring coil is the same for every design. Zips without one use `designs/Coil_shared.step`.

## Results

All results are per pole, on the L1 mesh, at 36 A unless stated. "3° window" means the average and minimum torque over the best 3° of stroke. "Per phase" means 40 × the per-pole value.

### Rotor tooth depth and back iron (`results/figures/compare_round3.png`)

| Design | Avg torque, 3° window (N·m) | Min torque, 3° window (N·m) | Peak torque (N·m) | Per phase (N·m) | Aligned/unaligned flux linkage | Iron per pole (cm³) | Torque per cm³ (N·m/cm³) |
|---|---|---|---|---|---|---|---|
| 10_10_5.6_10_3.5_2_3 (baseline) | 1.100 | 0.640 | 1.361 | 44.0 | 1.65 | 8.52 | 0.129 |
| **13_10_5.6_10_3.5_2_3** | **1.208** | 0.746 | 1.405 | 48.3 | 1.80 | 9.32 | 0.130 |
| 16_10_5.6_10_3.5_2_3 | 1.223 | 0.765 | 1.407 | 48.9 | 1.83 | 10.12 | 0.121 |
| 10_10_5.6_10_2.1_2_3 | 1.089 | 0.642 | 1.355 | 43.6 | 1.61 | 7.77 | 0.140 |
| 13_10_5.6_10_2.1_2_3 | 1.168 | 0.692 | 1.398 | 46.7 | 1.70 | 8.57 | 0.136 |

- **Rotor tooth depth.** Going from 10 to 13 mm adds 9.9% torque; 13 to 16 mm adds only 1.2% for 8.6% more iron. At aligned the 10 mm tooth runs at 1.91 T, the 13 mm at 1.72 T and the 16 mm at 1.46 T. Once the tooth stops saturating, the stator leg tips (about 1.9 T) limit the torque. The torque curve gets flatter on top from 10 to 13 mm. Above 3° all designs follow the same curve, which the tooth width sets.
- **Back iron.** Thinning it to 2.1 mm costs 1% with the 10 mm rotor tooth but 3.4% with the 13 mm one: the deeper tooth pushes more flux through the back iron, which then reaches 1.74 T.
- **Choice so far: 13_10_5.6_10_3.5_2_3.**

### Flux leakage at aligned (`results/figures/leakage_*.png`)

The share of the middle core's back-iron flux that passes through the middle rotor tooth at its mid-plane (cross-section integrals):

| Rotor tooth depth | Through the tooth | Leakage |
|---|---|---|
| 10 mm | 87.0% | 13.0% |
| 13 mm | 93.1% | 6.9% |
| 16 mm | 95.5% | 4.5% |

The figures draw only the field lines that do *not* go through the tooth. The largest remaining path is leg to leg across the coil window. With the 10 mm tooth there is also leakage around the axial ends, which mostly disappears with deeper teeth. Leakage to neighbouring poles is under 0.5%.

### Phase current, design 13_10_5.6_10_3.5_2_3 (`results/figures/current_sweep_13.png`)

Only 4 angles were run per current (0.5°, 1.5°, 2.5°, 3.5°), not full sweeps.

| Current | Current density | 0.5° | 1.5° | 2.5° | 3.5° | Mean of 4 angles vs 36 A | Copper loss vs 36 A |
|---|---|---|---|---|---|---|---|
| 36 A | 7.5 A/mm² | 0.733 | 1.284 | 1.405 | 0.797 | 1.00× | 1.00× |
| 45 A | 9.4 A/mm² | 0.843 | 1.636 | 1.963 | 1.205 | 1.34× | 1.56× |
| 50 A | 10.5 A/mm² | 0.891 | 1.790 | 2.223 | 1.432 | 1.50× | 1.93× |
| 60 A | 12.6 A/mm² | 0.970 | 2.033 | 2.624 | 1.829 | 1.77× | 2.78× |
| 70 A | 14.7 A/mm² | 1.032 | 2.223 | 2.914 | 2.133 | 1.97× | 3.78× |

Torque values are in N·m per pole.

- **Torque per amp stays about 30 mN·m/A from 36 to 70 A.** Torque per watt of copper loss halves over the same range.
- **Early angles saturate first.** At 0.5° the torque barely rises with current. So the flat top seen at 36 A becomes a peak at 2.5–3.5° at high current.
- **At 60–70 A the back iron (1.8–1.86 T) and the rotor tooth (1.9 T) saturate as well.** The best rotor depth and back-iron thickness are therefore likely thicker at high current than at 36 A.

## Caveats

- **B-H curve.** `scripts/m19_bh.dat` is a generic M-19 curve, not the sourced AK Steel / EMERF curve in `basic_mag_analysis/materials/M19_DiMax_BH.csv`. They agree within about 13% up to 1.3 T and within 7% above 1.8 T. In the 1.4–1.7 T knee the generic curve needs 20–30% less H, so absolute torque is probably a few percent optimistic. Comparisons between designs should hold. Switching to the sourced curve (extended above 2.1 T) is a planned fix.
- **Tooth-flux bug, fixed 2026-11-21.** `post_sector.py` sampled the rotor tooth only out to r = 112 mm, which under-counted `phi_tooth_uWb` for teeth deeper than 12 mm. The affected runs were reprocessed before these files were exported. Torque and the other quantities were never affected.
- **`MIN_FRAG`.** At 3° the segment edge leaves a thin sliver (about 2.5%) of the outermost rotor tooth. For the 16 mm tooth this made the linear solver diverge. `MIN_FRAG=0.05` drops slivers below 5% of a tooth. On a case that already converged it changed the torque by 0.04%. It is off by default and was used only for 16_10_5.6_10_3.5_2_3 at 3°.
- **Traced leakage figures overstate leakage.** Field lines are traced through B values averaged at the mesh nodes, so lines running close to the tooth surface drift out of it. The traced "useful" share comes out 8–10 points below the cross-section numbers above. Use the figures to see where leakage goes, and the table for how much.
- **L1 mesh only.** Results have not been re-checked on the finer L2 mesh. Before committing to a final design, re-run it at `LEVEL=L2`.

## Layout

```
designs/        input STEP zips, naming convention, shared ring coil
scripts/        sector.geo (geometry + mesh), case_sector.sif (Elmer case), m19_bh.dat,
                run_sector.sh (one angle), sweep_sector.sh (angle sweep), post_sector.py (per-run post),
                collect_sector.py / compare_designs.py (tables and plots), leakage_lines.py,
                probe.py / readvtu.py (VTU reading and point probing)
results/
  phase_<design>_L1[_A<amps>]/<angle>deg/summary.json   every quantity for one run
  phase_<design>_L1/results.csv, design_summary.json    per-design torque curve and summary
  figures/      comparison, leakage and current-sweep plots
  logs/         sweep logs
```

Meshes and field files (about 60 MB per run) are not committed. Re-running a case regenerates them under `runs/`.
