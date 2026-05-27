import argparse
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from bins import RAW_ARGS

ROOT = "/opt/conda/envs/root618/bin/root"


def arg_text(args):
    return ",".join(f"{x:g}" if isinstance(x, float) else str(x) for x in args)




def root_cmd(macro, args):
    return f"{ROOT} -l -b -q '{macro}({arg_text(args)})'"


def safe_log_name(stage, args):
    parts = [stage, "pt", args[0], args[1], "y", args[2], args[3], "cent", args[4], args[5]]
    if len(args) > 6:
        parts += ["cat", args[6]]
    name = "_".join(str(x).replace(".", "p") for x in parts)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name) + ".log"


def run_one(stage, macro, args, log_dir):
    cmd = root_cmd(macro, args)
    log_path = log_dir / safe_log_name(stage, args)
    with open(log_path, "w", encoding="utf-8", errors="ignore") as log:
        log.write(f"# STAGE: {stage}\n")
        log.write(f"# CMD: {cmd}\n\n")
        proc = subprocess.run(cmd, shell=True, stdout=log, stderr=subprocess.STDOUT)
    return {
        "stage": stage,
        "macro": macro,
        "args": args,
        "cmd": cmd,
        "returncode": proc.returncode,
        "log": log_path,
    }


def build_stages():
    raw_args = RAW_ARGS
    return [
        #("mc_massFit", [( "mc_MassFit_HighpT.C", args) for args in raw_args]),
        #("massFit", [( "MassFit_FixPar_Data.C", args) for args in raw_args]),
        ("ctauErr", [( "CtauErr.C", args) for args in raw_args]),
        ("ctauRes", [( "CtauRes.C", args) for args in raw_args]),
        ("ctauBkg", [( "CtauBkg_2exp.C", args) for args in raw_args]),
        ("ctauTrue", [( "CtauTrue.C", args) for args in raw_args]),
        ("final2D", [( "Final2DFit.C", args) for args in raw_args]),
    ]


def configure_root618_env():
    os.environ["CONDA_PREFIX"] = "/opt/conda/envs/root618"
    os.environ["PATH"] = f"/opt/conda/envs/root618/bin:{os.environ.get('PATH', '')}"
    os.environ["CPATH"] = (
        "/opt/conda/envs/root618/x86_64-conda-linux-gnu/sysroot/usr/include:"
        "/opt/conda/envs/root618/x86_64-conda-linux-gnu/sysroot/usr/include/x86_64-linux-gnu:"
        "/usr/include:/usr/include/x86_64-linux-gnu"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--workers", type=int, default=4)
    parser.add_argument("--log-dir", default="logs/pTreweight_all2D")
    args = parser.parse_args()

    configure_root618_env()
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    for stage, jobs in build_stages():
        print(f"\n===== {stage}: {len(jobs)} jobs, workers={args.workers} =====", flush=True)
        failures = []
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(run_one, stage, macro, job_args, log_dir) for macro, job_args in jobs]
            for future in as_completed(futures):
                result = future.result()
                status = "OK" if result["returncode"] == 0 else f"FAIL rc={result['returncode']}"
                print(f"[{status}] {result['macro']}({arg_text(result['args'])}) log={result['log']}", flush=True)
                if result["returncode"] != 0:
                    failures.append(result)

        if failures:
            print(f"\nSTOP: {stage} failed for {len(failures)} jobs")
            for failure in failures:
                print(f"- rc={failure['returncode']} log={failure['log']}\n  {failure['cmd']}")
            raise SystemExit(1)

    print("\nAll pTreweight 2D-fit stages finished.")


if __name__ == "__main__":
    main()
