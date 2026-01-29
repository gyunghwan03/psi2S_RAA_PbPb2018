#!/bin/bash

# 절대경로로 고정
BASE_DIR_ABS="$(cd "systematic_Pb" 2>/dev/null && pwd)"
if [[ -z "$BASE_DIR_ABS" ]]; then
  echo "❌ BASE_DIR systematic_Pb not found."
  exit 1
fi
LOG_ROOT="$BASE_DIR_ABS"   # 로그를 systematic_Pb 바로 아래에 저장
mkdir -p "$LOG_ROOT"

DIRS=(
  "HFup"
  "HFdown"
  "bkgPDF"
  "sigPDF"
  "sigPAR/SigPAR_alpha"
  "sigPAR/SigPAR_f"
  "sigPAR/SigPAR_n"
  "sigPAR/SigPAR_x"
)

start_time=$(date +%s)
LOG_PATHS=()

sanitize() { echo "${1//\//_}"; }

for DIR in "${DIRS[@]}"; do
  WORKDIR="${BASE_DIR_ABS}/${DIR}"

  echo "▶ Entering directory: $DIR"
  if [[ ! -d "$WORKDIR" ]]; then
    echo "⚠️  Skip: directory not found -> $WORKDIR"
    echo "--------------------------------------"
    continue
  fi
  cd "$WORKDIR" || { echo "❌ Failed to enter $WORKDIR"; echo "--------------------------------------"; continue; }

  base="$(sanitize "$DIR")"
  log_path="${LOG_ROOT}/${base}.log"   # 예: systematic_Pb/sigPAR_SigPar_f.log
  err_path="${LOG_ROOT}/${base}.err"
  tmp_path="${LOG_ROOT}/${base}.tmp"
  LOG_PATHS+=("$log_path")

  if [[ -x "./01mass.sh" ]]; then
    echo "🚀 Running 01mass.sh in $DIR..."
    sub_start=$(date +%s)

    # 전체 출력은 임시 파일(절대경로)에 수집
    ./01mass.sh > "$tmp_path" 2>&1
    rc=$?

    sub_end=$(date +%s)
    sub_elapsed=$(( sub_end - sub_start ))

    # 에러 키워드 탐지
    if [[ -f "$tmp_path" ]] && grep -qiE "error|segmentation violation|fatal" "$tmp_path"; then
      {
        echo "===== ERROR SUMMARY ($(date '+%F %T')) ====="
        grep -niE "error|segmentation violation|fatal" "$tmp_path" | tail -n 200
      } > "$err_path"
      mv -f "$tmp_path" "$log_path"   # 에러면 전체 로그 보존
      echo "⚠️  $DIR had errors (exit=$rc, ${sub_elapsed}s). Full log: $log_path, Summary: $err_path"
    else
      # 정상/무에러: 마지막 200줄만 보관(용량 절감)
      if [[ -f "$tmp_path" ]]; then
        tail -n 200 "$tmp_path" > "$log_path"
        rm -f "$tmp_path" "$err_path"
      else
        echo "⚠️  No tmp log created. (exit=$rc)"
        : > "$log_path"
      fi
      if (( rc == 0 )); then
        echo "✅ Finished $DIR successfully in ${sub_elapsed}s. (kept last 200 lines)"
      else
        echo "⚠️  $DIR exit=$rc but no error keywords. (kept last 200 lines)"
      fi
    fi
  else
    echo "⚠️  01mass.sh not found or not executable in $WORKDIR"
  fi

  cd - >/dev/null
  echo "--------------------------------------"
done

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))

echo ""
echo "🎯 All scripts completed. Total elapsed time: ${elapsed}s."
echo ""
echo "📄 Log files saved at:"
for path in "${LOG_PATHS[@]}"; do
  abs_path="$path"  # 이미 절대경로
  echo "   - $abs_path"
done
echo ""
