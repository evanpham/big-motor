#!/usr/bin/env bash
# Rotor-angle sweep of one design with the multi-pole phase model.
# usage: sweep_sector.sh <folder or .zip with the 3 STEP files> [angles...]   (default 0 0.5 ... 4.5)
# env: LEVEL=L1|L2 (default L1), NP (default 4), NSEC is not exposed (3 poles)
cd "$(dirname "$0")/.."
SRC=$1; shift
if [[ "$SRC" == *.zip ]]; then d="${SRC%.zip}"; [ -d "$d" ] || { mkdir "$d" && unzip -q "$SRC" -d "$d"; }; SRC=$d; fi
ANG=("$@"); [ ${#ANG[@]} -eq 0 ] && ANG=(0 0.5 1 1.5 2 2.5 3 3.5 4 4.5)
case "${LEVEL:=L1}" in
  L1) export H_GAP=0.5 H_IRON=1.2 H_COIL=2.5 ;;
  L2) export H_GAP=0.35 H_IRON=0.9 H_COIL=2.0 ;;
esac
export NP=${NP:-4}
name="phase_$(basename "$SRC")_$LEVEL"
for a in "${ANG[@]}"; do
  tag=$(printf "%04.1fdeg" "$a")
  [ -f "runs/$name/$tag/summary.json" ] && { echo "=== $tag already done"; continue; }
  echo "=== $tag  $(date +%H:%M:%S)"
  scripts/run_sector.sh "$SRC" "$a" "$name/$tag" | grep -E "WARNING|solve time|torque_Nm|failed" || echo "FAILED $a"
done
echo "=== sweep done $(date +%H:%M:%S)  -> runs/$name"
