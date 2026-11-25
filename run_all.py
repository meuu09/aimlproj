"""Rebuild every derived result from data/ in dependency order; dependents of a failed step are skipped.

Modes: default = everything offline except the original-code baseline; --baseline adds the baseline (first);
--quick = USGS-only analyses (no combined catalog, no robustness); --download also refreshes the downloads.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, "results", "run_log.txt")
ANALYSES = ["analysis.him_ides_himachal", "analysis.kishtwar_cluster", "analysis.depth_clusters"]
CASES = [[], ["--include-fixed"]]


def m(module, *args):
    """Command that runs a project script as a module from the repository root."""
    return ["-m", module, *args]


def steps(argv):
    """(name, command, dependencies) in dependency order for the chosen mode."""
    quick, download = "--quick" in argv, "--download" in argv
    s = []
    if "--baseline" in argv:
        s += [("baseline", m("baseline.reproduce_baseline", "--fig10"), []),
              ("baseline_comparison", m("baseline.baseline_comparison"), ["baseline"])]
    if download:
        s += [("usgs_depth_types", m("fetch.fetch_depth_type"), [])]
    if not quick:
        s += [("isc_parse", m("fetch.fetch_isc_bulletin", *([] if download else ["--parse-only"])), []),
              ("combine", m("catalog.build_combined_catalog"), ["isc_parse"] + (["usgs_depth_types"] if download else [])),
              ("qc", m("catalog.catalog_qc"), ["combine"]),
              ("completeness", m("robustness.mc_through_time"), ["combine"]),
              ("focal_mechanisms", m("fetch.fetch_focal_mechanisms", *([] if download else ["--parse-only"])),
               ["combine"])]
    cats = [([], "usgs", [])] if quick else [([], "usgs", []), (["--catalog", "combined"], "combined", ["combine"]),
                                             (["--catalog", "combined", "--drop-duplicates"], "nodup", ["combine"])]
    for cat_args, tag, deps in cats:
        for case in CASES:
            for module in ANALYSES:
                name = f"{module.split('.')[-1]}:{tag}{':incl' if case else ''}"
                s.append((name, m(module, *cat_args, *case), deps))
    if not quick:
        s += [("temporal", m("analysis.temporal_summary"), ["combine", "completeness"]),
              ("migration", m("analysis.migration_search"), ["combine"]),
              ("clustering", m("analysis.clustering_checks"), ["combine", "focal_mechanisms"]),
              ("bvalue_robustness", m("robustness.bvalue_robustness"), ["combine"]),
              ("cluster_robustness", m("robustness.cluster_robustness"), ["combine"]),
              ("depth_sensitivity", m("robustness.depth_sensitivity"), ["combine"])]
    s += [("tests", m("pytest", "-q", "tests"), [n for n, _, _ in s]),
          ("manifest", m("tools.write_manifest"), ["tests"])]
    return s


def run(args, log):
    """Run one step with the current Python; log its status and time."""
    t0 = time.time()
    r = subprocess.run([sys.executable] + args, cwd=HERE, capture_output=True, text=True, encoding="utf-8")
    status = "ok" if r.returncode == 0 else f"FAILED ({r.returncode})"
    line = f"{status:>12}  {time.time() - t0:6.1f}s  {' '.join(args)}"
    print(line, flush=True)
    log.write(line + "\n")
    if r.returncode:
        detail = (r.stdout[-2000:] + "\n" + r.stderr[-4000:]).strip()
        print(detail)
        log.write(detail + "\n")
    return r.returncode == 0


def main():
    """Run the steps; a step whose dependency failed or was skipped is skipped too."""
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    ok, bad = set(), []
    with open(LOG, "w", encoding="utf-8") as log:
        log.write(f"run_all.py {' '.join(sys.argv[1:])} started {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        for name, cmd, deps in steps(sys.argv[1:]):
            missing = [d for d in deps if d not in ok]
            if missing:
                line = f"{'SKIPPED':>12}  {'':>7}  {' '.join(cmd)}  (upstream failed: {', '.join(missing)})"
                print(line); log.write(line + "\n"); bad.append(name)
                continue
            (ok.add(name) if run(cmd, log) else bad.append(name))
        log.write(f"finished: {len(ok)} ok, {len(bad)} failed or skipped\n")
    print("all done" if not bad else f"{len(bad)} step(s) failed or skipped: {', '.join(bad)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
