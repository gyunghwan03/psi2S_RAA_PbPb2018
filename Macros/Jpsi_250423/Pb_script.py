import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

commands_nominal = [
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12, 15, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(15, 20, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(20, 25, 0, 1.6, 0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(25, 40, 0, 1.6, 0, 180, 1)'", 

"root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9, 12, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12, 15, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(15, 20, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(20, 25, 0, 1.6, 0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(25, 40, 0, 1.6, 0, 180, 2)'", 


"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 0, 20, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 20, 40, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 40, 60, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 60, 80, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 80, 100, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 100, 180, 1)'", 


"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 0, 20, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 20, 40, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 40, 60, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 60, 80, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 80, 100, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5, 40, 0, 1.6, 100, 180, 2)'", 


"root -l -b -q MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4,  0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5,   9, 1.6, 2.4,  0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9,    12, 1.6, 2.4,  0, 180, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12,   40, 1.6, 2.4,  0, 180, 1)'", 

"root -l -b -q MassFit_FixPar_Data.C'(3.5, 6.5, 1.6, 2.4,  0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(6.5,   9, 1.6, 2.4,  0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(9,    12, 1.6, 2.4,  0, 180, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(12,   40, 1.6, 2.4,  0, 180, 2)'", 

"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 0, 20, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 20, 60, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 60, 100, 1)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 100, 180, 1)'", 


"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 0, 20, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 20, 60, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 60, 100, 2)'", 
"root -l -b -q MassFit_FixPar_Data.C'(3.5, 40, 1.6, 2.4, 100, 180, 2)'", 
]

commands_pTreweight = [
    "root -l -b -q mc_MassFit_HighpT.C'(6.5, 7.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(7.5, 8.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(8.5, 9.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(9.5,  11, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(11,   13, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(13,   15, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(15, 17.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(17.5, 20, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(20, 22.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(22.5, 25, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(25, 27.5, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(27.5, 30, 0, 1.6, 0, 180, 1)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(30,   40, 0, 1.6, 0, 180, 1)'", 

    "root -l -b -q mc_MassFit_HighpT.C'(3.5, 4.5, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q mc_MassFit_HighpT.C'(4.5, 6.5, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q mc_MassFit_HighpT.C'(6.5, 8.5, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q mc_MassFit_HighpT.C'(8.5,  12, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q mc_MassFit_HighpT.C'(12,   15, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q mc_MassFit_HighpT.C'(15,   20, 1.6, 2.4, 0, 180, 1)'",
    "root -l -b -q mc_MassFit_HighpT.C'(20,   40, 1.6, 2.4, 0, 180, 1)'",

    "root -l -b -q mc_MassFit_HighpT.C'(6.5, 7.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(7.5, 8.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(8.5, 9.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(9.5,  11, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(11,   13, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(13,   15, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(15, 17.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(17.5, 20, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(20, 22.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(22.5, 25, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(25, 27.5, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(27.5, 30, 0, 1.6, 0, 180, 2)'", 
    "root -l -b -q mc_MassFit_HighpT.C'(30,   40, 0, 1.6, 0, 180, 2)'", 

    "root -l -b -q mc_MassFit_HighpT.C'(3.5, 4.5, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q mc_MassFit_HighpT.C'(4.5, 6.5, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q mc_MassFit_HighpT.C'(6.5, 8.5, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q mc_MassFit_HighpT.C'(8.5,  12, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q mc_MassFit_HighpT.C'(12,   15, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q mc_MassFit_HighpT.C'(15,   20, 1.6, 2.4, 0, 180, 2)'",
    "root -l -b -q mc_MassFit_HighpT.C'(20,   40, 1.6, 2.4, 0, 180, 2)'",
]



LOG_DIR = os.path.join(os.path.dirname(__file__), "logs")
os.makedirs(LOG_DIR, exist_ok=True)


def _parse_bin_label(cmd: str) -> str:
    """주어진 root 실행 커맨드 문자열에서 macro 이름과 (pt,y,cent,cat) 정보를 파싱해 사람이 읽을 수 있는 라벨로 반환.

    예시 입력:
      "root -l -b -q MassFit_FixPar_Data.C'(6.5, 9, 0, 1.6, 0, 180, 1)'"
    출력 라벨:
      "MassFit_FixPar_Data.C | pT[6.5,9], y[0,1.6], cent[0,180], cat=1"
    """
    try:
        # -q 다음 토큰부터 매크로와 인자를 포함
        q_idx = cmd.find("-q ")
        tail = cmd[q_idx + 3 :].strip() if q_idx != -1 else cmd
        # macro 이름: 따옴표 전까지
        macro = tail.split("'", 1)[0].strip()
        # 괄호 안 인자 추출
        lpar = tail.find("(")
        rpar = tail.rfind(")")
        args_str = tail[lpar + 1 : rpar]
        # 7개 인자 예상: pt1, pt2, y1, y2, cent1, cent2, cat
        parts = [p.strip() for p in args_str.split(",")]
        if len(parts) >= 7:
            pt1, pt2, y1, y2, c1, c2, cat = parts[:7]
            return f"{macro} | pT[{pt1},{pt2}], y[{y1},{y2}], cent[{c1},{c2}], cat={cat}"
        # 예외 시 원문 반환
        return tail
    except Exception:
        return cmd


def run_and_capture(cmd: str):
    """커맨드를 실행하고 결과(라벨, 리턴코드, [Stuck] 포함 여부, 로그 경로)를 반환."""
    label = _parse_bin_label(cmd)
    # 로그 파일명 안전화
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)[:200]
    log_path = os.path.join(LOG_DIR, f"{safe_name}.log")

    # 실행 및 출력 캡처(표준출력/표준에러 합침)
    start = datetime.now()
    try:
        proc = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        output = proc.stdout or ""
        rc = proc.returncode
    except Exception as e:
        output = f"[PythonException] {e}"
        rc = 1

    # 로그 저장
    try:
        with open(log_path, "w", encoding="utf-8", errors="ignore") as f:
            f.write(f"# CMD: {cmd}\n")
            f.write(f"# LABEL: {label}\n")
            f.write(f"# START: {start.isoformat()}\n")
            f.write(f"# END: {datetime.now().isoformat()}\n")
            f.write("\n")
            f.write(output)
    except Exception:
        # 로그 쓰기 실패는 무시
        pass

    # [Stuck] 탐지 (대소문자 구분 그대로, 추가로 소문자 stuck도 탐지)
    stuck = ("[Stuck]" in output) or ("[STUCK]" in output) or ("[stuck]" in output)

    return {
        "command": cmd,
        "label": label,
        "returncode": rc,
        "stuck": stuck,
        "log_path": log_path,
    }


def summarize(results):
    total = len(results)
    stucks = [r for r in results if r["stuck"]]
    failures = [r for r in results if r["returncode"] != 0]

    print("\n===== Fit summary =====")
    print(f"Total bins tried: {total}")
    print(f"[Stuck] detected: {len(stucks)}")
    print(f"Non-zero return code: {len(failures)}")

    if stucks:
        print("\n-- Bins with [Stuck] --")
        for r in stucks:
            print(f"- {r['label']}  (log: {os.path.relpath(r['log_path'], start=os.path.dirname(__file__))})")

    if failures:
        print("\n-- Bins with non-zero return code --")
        for r in failures:
            print(f"- rc={r['returncode']} | {r['label']}  (log: {os.path.relpath(r['log_path'], start=os.path.dirname(__file__))})")

    # CSV로도 저장
    csv_path = os.path.join(LOG_DIR, "fit_summary.csv")
    try:
        import csv

        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["label", "returncode", "stuck", "log_path"])            
            for r in results:
                w.writerow([r["label"], r["returncode"], int(r["stuck"]), r["log_path"]])
        print(f"\nSummary saved: {csv_path}")
    except Exception:
        pass


def main():
    # 여기서 어떤 커맨드 세트를 실행할지 선택
    commands = commands_pTreweight  # 필요 시 commands_nominal 등으로 교체/확장

    results = []
    print(f"Running {len(commands)} fits with up to 4 workers... Logs: {LOG_DIR}")
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_cmd = {executor.submit(run_and_capture, c): c for c in commands}
        for future in as_completed(future_to_cmd):
            res = future.result()
            results.append(res)
            # 짧은 진행 상황 출력
            status = "STUCK" if res["stuck"] else ("OK" if res["returncode"] == 0 else f"RC={res['returncode']}")
            print(f"[{status}] {res['label']}")

    summarize(results)


if __name__ == "__main__":
    main()
