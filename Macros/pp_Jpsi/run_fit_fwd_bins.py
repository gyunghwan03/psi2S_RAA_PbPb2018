#!/usr/bin/env python3
"""
Run the pp_Jpsi data fit chain for the 7 fwd y[1.6, 2.4] pT bins.

Chain (per bin):
  MassFit_FixPar_Data.C -> CtauErr.C -> CtauRes.C
  -> CtauBkg_LowPt.C -> CtauTrue.C -> Final2DFit.C

The MC mass fits (mc_MassFit_HighpT.C) are assumed to be already produced in
roots_MC/Mass/. CtauBkg_LowPt.C can be run in fast mode (FAST_CTAUBKG_NO_PLOT=1)
to skip the heavy plot/projection stage when only the fit result file is needed.

Defaults match the bin list used by compareDataToMC/make_ptw_aggregated_260519.C
for fwd rapidity, excluding pt5.0-6.5 (which already has a 2025-10-21 product).

Usage examples:
  python3 run_fit_fwd_bins.py                       # all 7 bins, full chain, sequential
  python3 run_fit_fwd_bins.py --bins 3.0-4.0 4.0-5.0
  python3 run_fit_fwd_bins.py --fast-bkg            # FAST_CTAUBKG_NO_PLOT=1
  python3 run_fit_fwd_bins.py --skip-existing       # skip a stage if its output exists
  python3 run_fit_fwd_bins.py --stages Mass CtauErr # subset of stages
  python3 run_fit_fwd_bins.py --jobs 2              # run two bins in parallel

All output is logged under logs/fit_fwd_pp_<UTC timestamp>/<bin>_<stage>.log
relative to the Macros/pp_Jpsi directory.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

DEFAULT_BINS: list[tuple[float, float]] = [
    (3.0, 4.0),
    (4.0, 5.0),
    (6.5, 8.5),
    (8.5, 12.0),
    (12.0, 15.0),
    (15.0, 20.0),
    (20.0, 40.0),
]

Y_LOW = 1.6
Y_HIGH = 2.4

DEFAULT_ROOT_BIN = "/opt/conda/envs/root618/bin/root"

# (stage_name, macro_filename, args_template, output_path_template)
# args are the ROOT macro arguments in (ptLow, ptHigh, yLow, yHigh, ...) order.
# Output path is relative to Macros/pp_Jpsi/. Used for --skip-existing.
STAGES: list[dict] = [
    {
        "name": "Mass",
        "macro": "MassFit_FixPar_Data.C",
        "args": "({ptL}, {ptH}, {yL}, {yH}, 1)",
        "out": "roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_pt{ptL_s}-{ptH_s}_y{yL_s}-{yH_s}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    },
    {
        "name": "CtauErr",
        "macro": "CtauErr.C",
        "args": "({ptL}, {ptH}, {yL}, {yH}, 1)",
        "out": "roots/2DFit_No_Weight/CtauErr/CtauErrResult_pt{ptL_s}-{ptH_s}_y{yL_s}-{yH_s}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    },
    {
        "name": "CtauRes",
        "macro": "CtauRes.C",
        "args": "({ptL}, {ptH}, {yL}, {yH}, 1)",
        "out": "roots/2DFit_No_Weight/CtauRes/CtauResResult_pt{ptL_s}-{ptH_s}_y{yL_s}-{yH_s}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    },
    {
        "name": "CtauBkg",
        "macro": "CtauBkg_LowPt.C",
        "args": "({ptL}, {ptH}, {yL}, {yH}, 1)",
        "out": "roots/2DFit_No_Weight/CtauBkg/CtauBkgResult_pt{ptL_s}-{ptH_s}_y{yL_s}-{yH_s}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    },
    {
        "name": "CtauTrue",
        "macro": "CtauTrue.C",
        "args": "({ptL}, {ptH}, {yL}, {yH}, 0.0, 2, 0.08)",
        "out": "roots/2DFit_No_Weight/CtauTrue/CtauTrueResult_Inclusive_pt{ptL_s}-{ptH_s}_y{yL_s}-{yH_s}_muPt0.0.root",
    },
    {
        "name": "Final",
        "macro": "Final2DFit.C",
        "args": "({ptL}, {ptH}, {yL}, {yH}, 1)",
        "out": "roots/2DFit_No_Weight/Final/2DFitResult_pt{ptL_s}-{ptH_s}_y{yL_s}-{yH_s}_muPt0.0_PRw_Effw0_Accw0_PtW0_TnP0.root",
    },
]

STAGE_BY_NAME = {s["name"]: s for s in STAGES}


def fmt_edge(x: float) -> str:
    # Match macro file naming: "%.1f" (e.g. 3.0, 12.0)
    return f"{x:.1f}"


@dataclass
class StageResult:
    bin_label: str
    stage: str
    returncode: int
    seconds: float
    skipped: bool
    log_path: Path


def stage_output(work_dir: Path, stage: dict, ptL: float, ptH: float) -> Path:
    return work_dir / stage["out"].format(
        ptL=ptL, ptH=ptH, yL=Y_LOW, yH=Y_HIGH,
        ptL_s=fmt_edge(ptL), ptH_s=fmt_edge(ptH),
        yL_s=fmt_edge(Y_LOW), yH_s=fmt_edge(Y_HIGH),
    )


def build_command(stage: dict, ptL: float, ptH: float, root_bin: str) -> list[str]:
    args = stage["args"].format(
        ptL=ptL, ptH=ptH, yL=Y_LOW, yH=Y_HIGH,
    )
    macro_call = f"{stage['macro']}{args}"
    return [root_bin, "-l", "-b", "-q", macro_call]


def run_stage(
    work_dir: Path,
    log_dir: Path,
    bin_pair: tuple[float, float],
    stage: dict,
    root_bin: str,
    fast_bkg: bool,
    skip_existing: bool,
) -> StageResult:
    ptL, ptH = bin_pair
    bin_label = f"pt{fmt_edge(ptL)}-{fmt_edge(ptH)}"
    log_path = log_dir / f"{bin_label}_{stage['name']}.log"

    out_path = stage_output(work_dir, stage, ptL, ptH)
    if skip_existing and out_path.exists():
        log_path.write_text(f"[skip] output exists: {out_path}\n")
        return StageResult(bin_label, stage["name"], 0, 0.0, True, log_path)

    cmd = build_command(stage, ptL, ptH, root_bin)
    env = os.environ.copy()
    if fast_bkg and stage["name"] == "CtauBkg":
        env["FAST_CTAUBKG_NO_PLOT"] = "1"

    t0 = time.monotonic()
    with log_path.open("w") as fh:
        fh.write(f"# cwd: {work_dir}\n")
        fh.write(f"# cmd: {' '.join(cmd)}\n")
        if fast_bkg and stage["name"] == "CtauBkg":
            fh.write("# env: FAST_CTAUBKG_NO_PLOT=1\n")
        fh.write(f"# started: {datetime.now(timezone.utc).isoformat()}\n")
        fh.flush()
        proc = subprocess.run(
            cmd,
            cwd=str(work_dir),
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
        )
    seconds = time.monotonic() - t0

    return StageResult(bin_label, stage["name"], proc.returncode, seconds, False, log_path)


def run_bin(
    work_dir: Path,
    log_dir: Path,
    bin_pair: tuple[float, float],
    stages: list[dict],
    root_bin: str,
    fast_bkg: bool,
    skip_existing: bool,
    stop_on_fail: bool,
) -> list[StageResult]:
    results: list[StageResult] = []
    for stage in stages:
        res = run_stage(work_dir, log_dir, bin_pair, stage, root_bin, fast_bkg, skip_existing)
        tag = "skip" if res.skipped else ("ok " if res.returncode == 0 else "FAIL")
        print(
            f"  [{res.bin_label} {res.stage:8s}] {tag} "
            f"rc={res.returncode} time={res.seconds:6.1f}s log={res.log_path}",
            flush=True,
        )
        results.append(res)
        if stop_on_fail and res.returncode != 0 and not res.skipped:
            print(
                f"  [{res.bin_label}] stopping further stages: {res.stage} failed",
                flush=True,
            )
            break
    return results


def parse_bin(s: str) -> tuple[float, float]:
    try:
        lo_s, hi_s = s.split("-")
        return float(lo_s), float(hi_s)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"bin '{s}' must look like 'LOW-HIGH', e.g. 3.0-4.0"
        ) from exc


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--bins",
        nargs="+",
        type=parse_bin,
        default=DEFAULT_BINS,
        metavar="LOW-HIGH",
        help="pT bins (e.g. 3.0-4.0 4.0-5.0). Default: all 7 fwd bins.",
    )
    p.add_argument(
        "--stages",
        nargs="+",
        choices=[s["name"] for s in STAGES],
        default=[s["name"] for s in STAGES],
        help="Stages to run, in order. Default: full chain.",
    )
    p.add_argument(
        "--root-bin",
        default=DEFAULT_ROOT_BIN,
        help=f"Path to root executable. Default: {DEFAULT_ROOT_BIN}",
    )
    p.add_argument(
        "--fast-bkg",
        action="store_true",
        help="Set FAST_CTAUBKG_NO_PLOT=1 for the CtauBkg stage (skips plots).",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip a stage if its output ROOT already exists.",
    )
    p.add_argument(
        "--no-stop-on-fail",
        dest="stop_on_fail",
        action="store_false",
        help="Continue subsequent stages even if a stage fails (default: stop).",
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of bins to run in parallel. Each bin runs its stages serially. Default: 1.",
    )
    p.add_argument(
        "--work-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Working directory (must contain the .C macros). Default: directory of this script.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    work_dir: Path = args.work_dir.resolve()
    if not (work_dir / "MassFit_FixPar_Data.C").exists():
        print(f"ERROR: macros not found in {work_dir}", file=sys.stderr)
        return 2
    if not Path(args.root_bin).exists():
        print(f"ERROR: root binary not found: {args.root_bin}", file=sys.stderr)
        return 2

    stages = [STAGE_BY_NAME[name] for name in args.stages]

    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_dir = work_dir / "logs" / f"fit_fwd_pp_{ts}"
    log_dir.mkdir(parents=True, exist_ok=True)

    print(f"work_dir   : {work_dir}")
    print(f"root       : {args.root_bin}")
    print(f"stages     : {[s['name'] for s in stages]}")
    print(f"bins       : {[(fmt_edge(a), fmt_edge(b)) for a, b in args.bins]}")
    print(f"fast_bkg   : {args.fast_bkg}")
    print(f"skip_exist : {args.skip_existing}")
    print(f"jobs       : {args.jobs}")
    print(f"logs       : {log_dir}")
    print()

    all_results: list[StageResult] = []
    if args.jobs <= 1:
        for bin_pair in args.bins:
            print(f"=== bin pt{fmt_edge(bin_pair[0])}-{fmt_edge(bin_pair[1])}_y1.6-2.4 ===", flush=True)
            all_results.extend(
                run_bin(
                    work_dir, log_dir, bin_pair, stages,
                    args.root_bin, args.fast_bkg, args.skip_existing, args.stop_on_fail,
                )
            )
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            fut_to_bin = {
                ex.submit(
                    run_bin,
                    work_dir, log_dir, bin_pair, stages,
                    args.root_bin, args.fast_bkg, args.skip_existing, args.stop_on_fail,
                ): bin_pair
                for bin_pair in args.bins
            }
            for fut in as_completed(fut_to_bin):
                bp = fut_to_bin[fut]
                print(f"=== finished bin pt{fmt_edge(bp[0])}-{fmt_edge(bp[1])}_y1.6-2.4 ===", flush=True)
                all_results.extend(fut.result())

    # Summary
    print()
    print("=== Summary ===")
    fail = 0
    for r in all_results:
        tag = "SKIP" if r.skipped else ("OK  " if r.returncode == 0 else "FAIL")
        if not r.skipped and r.returncode != 0:
            fail += 1
        print(f"  {tag}  {r.bin_label:14s} {r.stage:8s} rc={r.returncode} t={r.seconds:6.1f}s")
    print(f"failed stages: {fail}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
