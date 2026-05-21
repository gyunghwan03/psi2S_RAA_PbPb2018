#!/usr/bin/env bash
# Full-statistics Eff/Acc for the H_NP_rational PtW candidate
# (NP: rational fit; PR: 2-exp symlinked from G_aggregateHighPt_2exp).
#
# Uses the Eff_Acc_260519_v2 macros and passes the ptw dir explicitly.
# Outputs land under Eff_Acc_260519_v2/roots/H_NP_rational_full/.
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
# Reuse 260515 macros (already compiled / parameterized).
MACRO_DIR="${SCRIPT_DIR}"
PTW_DIR="${REPO_DIR}/compareDataToMC/ptw_candidates_260520_v3/H_NP_rational"
ROOTS_OUT="${SCRIPT_DIR}/roots"
mkdir -p "${ROOTS_OUT}"

# Macros write ROOTs to <macro_dir>/roots. We will mv them after
# each run into our 260519_v2 roots/H_NP_rational_full directory.
TARGET_DIR="${ROOTS_OUT}/H_NP_rational_full"
mkdir -p "${TARGET_DIR}"

LOG_DIR="${REPO_DIR}/logs_summary/260519_v2_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_DIR}"
echo "[INFO] logs: ${LOG_DIR}"

fail=0
TAG="260520v3_H_NP_rational"

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

run_one "G_eff_pbpb_pr" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(-1,true,true,true,true,true,true,false,false,\"${TAG}\",\"${PTW_DIR}\")"
run_one "G_eff_pbpb_np" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(-1,false,true,true,true,true,true,false,false,\"${TAG}\",\"${PTW_DIR}\")"
run_one "G_eff_pp_pr" \
  "get_Eff_Jpsi_pp_hwan_v6.C+(-1,true,true,true,true,false,true,\"${TAG}\",\"${PTW_DIR}\")"
run_one "G_eff_pp_np" \
  "get_Eff_Jpsi_pp_hwan_v6.C+(-1,false,true,true,true,false,true,\"${TAG}\",\"${PTW_DIR}\")"
run_one "G_acc_pp_pr" \
  "JpsiaccStudy_v2.C+(-1,true,false,true,true,\"${TAG}\",\"${PTW_DIR}\")"
run_one "G_acc_pp_np" \
  "JpsiaccStudy_v2.C+(-1,false,false,true,true,\"${TAG}\",\"${PTW_DIR}\")"
run_one "G_acc_pbpb_pr" \
  "JpsiaccStudy_v2.C+(-1,true,true,true,true,\"${TAG}\",\"${PTW_DIR}\")"
run_one "G_acc_pbpb_np" \
  "JpsiaccStudy_v2.C+(-1,false,true,true,true,\"${TAG}\",\"${PTW_DIR}\")"

# Move outputs to our 260519_v2 dir.
for f in "${MACRO_DIR}/roots/"*"${TAG}"*.root; do
  [[ -e "${f}" ]] || continue
  mv -f "${f}" "${TARGET_DIR}/"
done
echo "[INFO] staged -> ${TARGET_DIR}"

if [[ "${fail}" -eq 0 ]]; then
  echo "[DONE] all 8 stages OK"
else
  echo "[DONE] ${fail} stage(s) failed"
  exit 1
fi
