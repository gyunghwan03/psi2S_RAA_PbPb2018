#!/usr/bin/env bash
# Run v6/v2 full-stat Eff/Acc with the H_NP_rational PtW candidate
# (option-2 aggregated bins + rational fit for NP, PR same as G_aggregate).
# Only PtW-on phase is needed (noPtW is candidate-independent).
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MACRO_DIR="${REPO_DIR}/Eff_Acc_260513"
PTW_DIR="${REPO_DIR}/compareDataToMC/ptw_candidates_260520_v3/H_NP_rational"
ROOTS_OUT="${SCRIPT_DIR}/roots"
mkdir -p "${ROOTS_OUT}"

LOG_DIR="${REPO_DIR}/logs_summary/260520_H_NP_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_DIR}"
echo "[INFO] logs: ${LOG_DIR}"
echo "[INFO] ptw dir: ${PTW_DIR}"

fail=0
TAG="260520_H_NP_rational"

run_one() {
  local label="$1" call="$2"
  local log="${LOG_DIR}/${label}.log"
  echo "[RUN] ${label}"
  if (cd "${MACRO_DIR}" && root -l -b -q "${call}") >"${log}" 2>&1; then
    echo "[OK ] ${label}"
  else
    echo "[ERR] ${label}  (log: ${log})"
    fail=$((fail+1))
  fi
}

# pp Eff (PR / NP) PtW on
run_one "H_eff_pp_pr" \
  "get_Eff_Jpsi_pp_hwan_v6.C+(1, true, 0, 6.5, 50.0, 0.0, 2.4, \"${PTW_DIR}\", \"${TAG}\")"
run_one "H_eff_pp_np" \
  "get_Eff_Jpsi_pp_hwan_v6.C+(2, true, 0, 6.5, 50.0, 0.0, 2.4, \"${PTW_DIR}\", \"${TAG}\")"
# PbPb Eff (PR / NP) PtW on
run_one "H_eff_pbpb_pr" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(1, true, 0, true, 3.0, 50.0, 0.0, 2.4, 0, 180, -1L, true, true, true, false, false, \"${PTW_DIR}\", \"${TAG}\")"
run_one "H_eff_pbpb_np" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(2, true, 0, true, 3.0, 50.0, 0.0, 2.4, 0, 180, -1L, true, true, true, false, false, \"${PTW_DIR}\", \"${TAG}\")"
# Acceptance (4 jobs) PtW on
run_one "H_acc_pp_pr"   "JpsiaccStudy_v2.C+(1, 1, 0, TString(\"${TAG}\"), -1L, false, true, TString(\"${PTW_DIR}\"))"
run_one "H_acc_pp_np"   "JpsiaccStudy_v2.C+(2, 1, 0, TString(\"${TAG}\"), -1L, false, true, TString(\"${PTW_DIR}\"))"
run_one "H_acc_pbpb_pr" "JpsiaccStudy_v2.C+(1, 1, 0, TString(\"${TAG}\"), -1L, true,  true, TString(\"${PTW_DIR}\"))"
run_one "H_acc_pbpb_np" "JpsiaccStudy_v2.C+(2, 1, 0, TString(\"${TAG}\"), -1L, true,  true, TString(\"${PTW_DIR}\"))"

DEST="${ROOTS_OUT}/H_NP_rational_full"
mkdir -p "${DEST}"
for f in "${MACRO_DIR}/roots/"*"${TAG}"*.root; do
  [[ -e "${f}" ]] || continue
  mv -f "${f}" "${DEST}/"
done
echo "[INFO] staged -> ${DEST}"

if [[ "${fail}" -eq 0 ]]; then
  echo "[DONE] all 8 stages OK"
else
  echo "[DONE] ${fail} stage(s) failed"
  exit 1
fi
