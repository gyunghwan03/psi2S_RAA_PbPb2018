#!/bin/bash

BASE_DIR="systematic_pp"
LOG_ROOT="${BASE_DIR}"   # logs/ 서브폴더 대신 BASE_DIR 바로 아래에 저장
mkdir -p "$LOG_ROOT"

DIRS=(
  "bkgPDF"
  "sigPDF"
  "sigPAR/SigPar_alpha"
  "sigPAR/SigPar_f"
  "sigPAR/SigPar_n"
  "sigPAR/SigPar_x"
)

start_time=$(date +%s)
LOG_PATHS=()

sanitize() { echo "${1//\//_}"; }

for DIR in "${DIRS[@]}"; do
  echo "▶ Entering directory: $DIR"
  cd "$BASE_DIR/$DIR" || { echo "❌ Failed to enter $BASE_DIR/$DIR"; exit 1; }

  base="$(sanitize "$DIR")"
  log_path="${LOG_ROOT}/${base}.log"   # 예: systematic_pp/sigPAR_SigPar_f.log
  err_path="${LOG_ROOT}/${base}.err"
  tmp_path="${LOG_ROOT}/${base}.tmp"
  LOG_PATHS+=("$log_path")

  if [ -x "./01mass.sh" ]; then
    echo "🚀 Running 01mass.sh in $DIR..."
    sub_start=$(date +%s)

    # 전체 출력은 임시 파일에 수집
    ./01mass.sh > "$tmp_path" 2>&1
    rc=$?

    sub_end=$(date +%s)
    sub_elapsed=$(( sub_end - sub_start ))

    # 에러 키워드 탐지
    if grep -qiE "error|segmentation violation|fatal" "$tmp_path"; then
      {
        echo "===== ERROR SUMMARY ($(date '+%F %T')) ====="
        grep -niE "error|segmentation violation|fatal" "$tmp_path" | tail -n 200
      } > "$err_path"
      mv -f "$tmp_path" "$log_path"   # 에러면 전체 로그 보존
      echo "⚠️  $DIR had errors (exit=$rc, ${sub_elapsed}s). Full log: $log_path, Summary: $err_path"
    else
      # 정상/무에러: 마지막 200줄만 보관(용량 절감)
      tail -n 200 "$tmp_path" > "$log_path"
      rm -f "$tmp_path" "$err_path"
      if [ $rc -eq 0 ]; then
        echo "✅ Finished $DIR successfully in ${sub_elapsed}s. (kept last 200 lines)"
      else
        echo "⚠️  $DIR exit=$rc but no error keywords. (kept last 200 lines)"
      fi
    fi
  else
    echo "⚠️  01mass.sh not found or not executable in $DIR"
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
  abs_path="$(realpath "$path" 2>/dev/null || echo "$path")"
  echo "   - $abs_path"
done
echo ""