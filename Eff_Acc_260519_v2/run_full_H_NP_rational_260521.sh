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
MACRO_DIR="${SCRIPT_DIR}"
ROOT_EXE=/opt/conda/envs/root634/bin/root
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

declare -A job_pids
declare -a job_labels

start_one() {
  local label="$1" call="$2"
  local log="${LOG_DIR}/${label}.log"
  echo "[RUN] ${label}"
  (cd "${MACRO_DIR}" && ${ROOT_EXE} -l -b -q "${call}") >"${log}" 2>&1 &
  job_pids["${label}"]=$!
  job_labels+=("${label}")
}

start_one "H_eff_pbpb_pr" \
  "get_Eff_JPsi_pbpb_hwan_v6.C(1,true,0,true,3.0,50.0,0.0,2.4,0,180,-1L,true,true,true,false,false,\"${PTW_DIR}\",\"${TAG}\")"
start_one "H_eff_pbpb_np" \
  "get_Eff_JPsi_pbpb_hwan_v6.C(2,true,0,true,3.0,50.0,0.0,2.4,0,180,-1L,true,true,true,false,false,\"${PTW_DIR}\",\"${TAG}\")"
start_one "H_eff_pp_pr" \
  "get_Eff_Jpsi_pp_hwan_v6.C(1,true,0,6.5,50.0,0.0,2.4,\"${PTW_DIR}\",\"${TAG}\")"
start_one "H_eff_pp_np" \
  "get_Eff_Jpsi_pp_hwan_v6.C(2,true,0,6.5,50.0,0.0,2.4,\"${PTW_DIR}\",\"${TAG}\")"
start_one "H_acc_pp_pr" \
  "JpsiaccStudy_v2.C(1,1,0,TString(\"${TAG}\"),-1L,false,true,TString(\"${PTW_DIR}\"))"
start_one "H_acc_pp_np" \
  "JpsiaccStudy_v2.C(2,1,0,TString(\"${TAG}\"),-1L,false,true,TString(\"${PTW_DIR}\"))"
start_one "H_acc_pbpb_pr" \
  "JpsiaccStudy_v2.C(1,1,0,TString(\"${TAG}\"),-1L,true,true,TString(\"${PTW_DIR}\"))"
start_one "H_acc_pbpb_np" \
  "JpsiaccStudy_v2.C(2,1,0,TString(\"${TAG}\"),-1L,true,true,TString(\"${PTW_DIR}\"))"

echo "[INFO] waiting for all 8 jobs..."
for label in "${job_labels[@]}"; do
  if wait "${job_pids[${label}]}"; then
    echo "[OK ] ${label}"
  else
    echo "[ERR] ${label}  (log: ${LOG_DIR}/${label}.log)"
    fail=$((fail+1))
  fi
done

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
