#!/usr/bin/env bash
# Full-statistics Eff/Acc using the legacy Eff_Acc_260513 v6/v2 macros:
#   - get_Eff_Jpsi_pp_hwan_v6.C        (pp efficiency)
#   - get_Eff_JPsi_pbpb_hwan_v6.C      (PbPb efficiency)
#   - JpsiaccStudy_v2.C                (acceptance)
#
# Two phases per call: PtW on (G_aggregateHighPt_2exp) and PtW off.
# Outputs land under Eff_Acc_260519_v2/roots/{G_aggregate_full,noPtW_full}.
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MACRO_DIR="${REPO_DIR}/Eff_Acc_260513"
PTW_DIR="${REPO_DIR}/compareDataToMC/ptw_candidates_260519_v2/G_aggregateHighPt_2exp"
ROOTS_OUT="${SCRIPT_DIR}/roots"
mkdir -p "${ROOTS_OUT}"

LOG_DIR="${REPO_DIR}/logs_summary/260520_v6v2_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_DIR}"
echo "[INFO] logs: ${LOG_DIR}"
echo "[INFO] ptw dir: ${PTW_DIR}"

fail=0

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

# Run a function defined inside an already-existing .C macro by emitting an
# inline driver file and pointing root at it. Used for Run_Eff_Jpsi_pp_hwan_v6
# (the internal "no outer wrapper" entry).
run_inline() {
  local label="$1" loadFile="$2" callExpr="$3"
  local log="${LOG_DIR}/${label}.log"
  local drv
  drv="$(mktemp /tmp/run_${label}_XXXXXX.C)"
  cat > "${drv}" <<EOF
{
  gROOT->LoadMacro("${loadFile}+");
  ${callExpr};
}
EOF
  echo "[RUN] ${label}"
  if (cd "${MACRO_DIR}" && root -l -b -q "${drv}") >"${log}" 2>&1; then
    echo "[OK ] ${label}"
  else
    echo "[ERR] ${label}  (log: ${log})"
    fail=$((fail+1))
  fi
  rm -f "${drv}"
}

stage_outputs() {
  local tag="$1" dest="$2"
  mkdir -p "${dest}"
  for f in "${MACRO_DIR}/roots/"*"${tag}"*.root; do
    [[ -e "${f}" ]] || continue
    mv -f "${f}" "${dest}/"
  done
  echo "[INFO] staged ${tag} -> ${dest}"
}

# =========================================================
# Phase 1 : G_aggregate full-stat, PtW ON
# =========================================================
echo
echo "===== Phase 1: PtW ON  (G_aggregateHighPt_2exp) ====="
TAG_ON="260519v2_G_aggregate"

# pp Eff (PR / NP) PtW on : outer wrapper get_Eff_Jpsi_pp_hwan_v6 always
# enables PtW, so we use that directly.
run_one "v6_eff_pp_pr" \
  "get_Eff_Jpsi_pp_hwan_v6.C+(1, true, 0, 6.5, 50.0, 0.0, 2.4, \"${PTW_DIR}\", \"${TAG_ON}\")"
run_one "v6_eff_pp_np" \
  "get_Eff_Jpsi_pp_hwan_v6.C+(2, true, 0, 6.5, 50.0, 0.0, 2.4, \"${PTW_DIR}\", \"${TAG_ON}\")"

# PbPb Eff (PR / NP) PtW on (isPtWeight=0)
run_one "v6_eff_pbpb_pr" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(1, true, 0, true, 3.0, 50.0, 0.0, 2.4, 0, 180, -1L, true, true, true, false, false, \"${PTW_DIR}\", \"${TAG_ON}\")"
run_one "v6_eff_pbpb_np" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(2, true, 0, true, 3.0, 50.0, 0.0, 2.4, 0, 180, -1L, true, true, true, false, false, \"${PTW_DIR}\", \"${TAG_ON}\")"

# Acceptance (4 jobs) PtW on (wtopt=1)
run_one "v2_acc_pp_pr"   "JpsiaccStudy_v2.C+(1, 1, 0, TString(\"${TAG_ON}\"), -1L, false, true, TString(\"${PTW_DIR}\"))"
run_one "v2_acc_pp_np"   "JpsiaccStudy_v2.C+(2, 1, 0, TString(\"${TAG_ON}\"), -1L, false, true, TString(\"${PTW_DIR}\"))"
run_one "v2_acc_pbpb_pr" "JpsiaccStudy_v2.C+(1, 1, 0, TString(\"${TAG_ON}\"), -1L, true,  true, TString(\"${PTW_DIR}\"))"
run_one "v2_acc_pbpb_np" "JpsiaccStudy_v2.C+(2, 1, 0, TString(\"${TAG_ON}\"), -1L, true,  true, TString(\"${PTW_DIR}\"))"

stage_outputs "${TAG_ON}" "${ROOTS_OUT}/G_aggregate_full"

# =========================================================
# Phase 2 : Full-stat, PtW OFF
# =========================================================
echo
echo "===== Phase 2: PtW OFF (noPtW) ====="
TAG_OFF="260518full_noPtW"

# pp Eff PtW off : need isPtW=false which the outer wrapper hard-codes to
# true. So we call the internal Run_Eff_Jpsi_pp_hwan_v6 directly via an
# inline driver script.
run_inline "v6_noPtW_eff_pp_pr" "get_Eff_Jpsi_pp_hwan_v6.C" \
  "Run_Eff_Jpsi_pp_hwan_v6(-1L, true,  true, true, false, true, \"\", \"${TAG_OFF}\")"
run_inline "v6_noPtW_eff_pp_np" "get_Eff_Jpsi_pp_hwan_v6.C" \
  "Run_Eff_Jpsi_pp_hwan_v6(-1L, false, true, true, false, true, \"\", \"${TAG_OFF}\")"

# PbPb Eff PtW off (isPtWeight = -999)
run_one "v6_noPtW_eff_pbpb_pr" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(1, true, -999, true, 3.0, 50.0, 0.0, 2.4, 0, 180, -1L, true, true, true, false, false, \"\", \"${TAG_OFF}\")"
run_one "v6_noPtW_eff_pbpb_np" \
  "get_Eff_JPsi_pbpb_hwan_v6.C+(2, true, -999, true, 3.0, 50.0, 0.0, 2.4, 0, 180, -1L, true, true, true, false, false, \"\", \"${TAG_OFF}\")"

# Acceptance PtW off (wtopt=0)
run_one "v2_noPtW_acc_pp_pr"   "JpsiaccStudy_v2.C+(1, 0, 0, TString(\"${TAG_OFF}\"), -1L, false, true, TString(\"\"))"
run_one "v2_noPtW_acc_pp_np"   "JpsiaccStudy_v2.C+(2, 0, 0, TString(\"${TAG_OFF}\"), -1L, false, true, TString(\"\"))"
run_one "v2_noPtW_acc_pbpb_pr" "JpsiaccStudy_v2.C+(1, 0, 0, TString(\"${TAG_OFF}\"), -1L, true,  true, TString(\"\"))"
run_one "v2_noPtW_acc_pbpb_np" "JpsiaccStudy_v2.C+(2, 0, 0, TString(\"${TAG_OFF}\"), -1L, true,  true, TString(\"\"))"

stage_outputs "${TAG_OFF}" "${ROOTS_OUT}/noPtW_full"

echo
if [[ "${fail}" -eq 0 ]]; then
  echo "[DONE] all 16 stages OK"
else
  echo "[DONE] ${fail} stage(s) failed"
  exit 1
fi
