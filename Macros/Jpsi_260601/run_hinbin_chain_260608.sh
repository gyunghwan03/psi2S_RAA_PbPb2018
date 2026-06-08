#!/usr/bin/env bash
# Full 2D-fit chain for HIN low-pT bins (PbPb fwd y1.6-2.4 cent0-180).
# Usage: run_hinbin_chain_260608.sh "3.0,4.5" [more bins...]
set -u
ROOT_EXE=/opt/conda/envs/root634/bin/root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
LOG="${SCRIPT_DIR}/hinbin_chain_260608.log"
: > "${LOG}"

run() {
  local label="$1" call="$2"
  echo "[$(date +%H:%M:%S)] START ${label}: ${call}" | tee -a "${LOG}"
  local t0=$SECONDS
  if ${ROOT_EXE} -l -b -q "${call}" >>"${LOG}" 2>&1; then
    echo "[$(date +%H:%M:%S)] OK    ${label}  ($((SECONDS-t0))s)" | tee -a "${LOG}"
  else
    echo "[$(date +%H:%M:%S)] FAIL  ${label}" | tee -a "${LOG}"
    return 1
  fi
}

for bin in "$@"; do
  IFS=',' read lo hi <<< "${bin}"
  tag="pt${lo}-${hi}"
  echo "===== BIN ${tag} =====" | tee -a "${LOG}"
  run "${tag} mcMass"  "mc_MassFit_HighpT.C(${lo},${hi},1.6,2.4,0,180)"   || exit 1
  run "${tag} datMass" "MassFit_FixPar_Data.C(${lo},${hi},1.6,2.4,0,180)" || exit 1
  run "${tag} ctauErr" "CtauErr.C(${lo},${hi},1.6,2.4,0,180)"             || exit 1
  run "${tag} ctauRes" "CtauRes.C(${lo},${hi},1.6,2.4,0,180)"             || exit 1
  run "${tag} ctauBkg" "CtauBkg_LowPt.C(${lo},${hi},1.6,2.4,0,180)"       || exit 1
  run "${tag} ctauTru" "CtauTrue.C(${lo},${hi},1.6,2.4,0,180)"            || exit 1
  run "${tag} final"   "Final2DFit.C(${lo},${hi},1.6,2.4,0,180)"          || exit 1
done
echo "[$(date +%H:%M:%S)] ALL DONE" | tee -a "${LOG}"
