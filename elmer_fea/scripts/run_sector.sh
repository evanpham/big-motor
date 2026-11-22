#!/usr/bin/env bash
# One rotor angle of the periodic-sector model.
# usage: run_sector.sh <folder with the 3 STEP files> <theta deg> <run name>
# env: H_GAP H_IRON H_COIL H_AIR (mm), NP (default 1), AMPS, TURNS
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC=$(cd "$1" && pwd); TH=$2; RUN="$ROOT/runs/$3"; NP=${NP:-1}
core=$(ls "$SRC"/*Stator_C_Core*.step); bar=$(ls "$SRC"/*Rotor_I_Bar*.step)
# the ring coil is the same for every design; use the shared one when a design zip has none
coil=$(ls "$SRC"/*Coil*.step 2>/dev/null || ls "$ROOT"/designs/Coil_shared.step)
mkdir -p "$RUN"; cd "$RUN"; rm -rf results mesh
G=(-setnumber theta "$TH")
[ -n "${NSEC:-}" ] && G+=(-setnumber nsec "$NSEC")
[ -n "${MIN_FRAG:-}" ] && G+=(-setnumber min_frag "$MIN_FRAG")
for v in GAP IRON COIL AIR; do n=H_$v; [ -n "${!n:-}" ] && G+=(-setnumber "h_${v,,}" "${!n}e-3"); done
gmsh -3 -nt 4 "${G[@]}" -setstring core "$core" -setstring coil "$coil" -setstring bar "$bar" -o mesh.msh "$ROOT/scripts/sector.geo" > gmsh.log 2>&1 || { grep Error gmsh.log | head; exit 1; }
grep -E "WARNING|volumes:|exterior" gmsh.log | sed 's/Info    : //'
A=$(sed -n 's/.*coil_area=//p' gmsh.log)
echo "  mesh: $(grep -E 'nodes .* elements' gmsh.log | tail -1 | sed 's/Info    : //'), coil area $(python3 -c "print(round($A*1e6,2))") mm2"
if [ "$NP" -gt 1 ]; then
  ElmerGrid 14 2 mesh.msh -autoclean -out mesh -metis "$NP" 4 > elmergrid.log 2>&1
else
  ElmerGrid 14 2 mesh.msh -autoclean -out mesh > elmergrid.log 2>&1
fi
cp "$ROOT/scripts/case_sector.sif" case.sif; cp "$ROOT/scripts/m19_bh.dat" .
sed -i "s/^\$ coilA = [^ ]*/\$ coilA = $A/" case.sif
[ -n "${AMPS:-}" ]  && sed -i "s/^\$ amps  = .*/\$ amps  = $AMPS/" case.sif
[ -n "${TURNS:-}" ] && sed -i "s/^\$ turns = .*/\$ turns = $TURNS/" case.sif
SECONDS=0
if [ "$NP" -gt 1 ]; then mpirun -np "$NP" ElmerSolver_mpi case.sif > solver.log 2>&1 || true
else ElmerSolver case.sif > solver.log 2>&1 || true; fi
grep -q "ALL DONE" solver.log || { echo "solver failed, see $RUN/solver.log"; exit 1; }
echo "  solve time: ${SECONDS}s on $NP process(es), nonlinear iterations: $(grep -c 'NS (ITER=.*mgdynamics$' solver.log)"
python3 -I "$ROOT/scripts/post_sector.py" "$RUN" "$TH" "$(sed -n 's/^\$ turns = //p' case.sif)" "$(sed -n 's/^\$ amps  = //p' case.sif)" "$A"
python3 -I "$ROOT/scripts/to_paraview.py" "$RUN/results" "$RUN/results/fields_$(basename "$3").vtu" > /dev/null
