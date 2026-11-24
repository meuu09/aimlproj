"""Rebuild every derived result from data/ (--quick: USGS only, --baseline: also the original code)."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["him_ides_himachal.py", "kishtwar_cluster.py", "depth_clusters.py"]
CATALOGS = [[], ["--catalog", "combined"], ["--catalog", "combined", "--drop-duplicates"]]
if "--quick" in sys.argv:
    CATALOGS = CATALOGS[:1]


def run(args):
    """Run one script with the current Python and report status and time."""
    t0 = time.time()
    r = subprocess.run([sys.executable] + args, cwd=HERE, capture_output=True, text=True, encoding="utf-8")
    status = "ok" if r.returncode == 0 else f"FAILED ({r.returncode})"
    print(f"{status:>10}  {time.time() - t0:5.1f}s  {' '.join(args)}", flush=True)
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-4000:], sep="\n")
    return r.returncode


def main():
    """Run every step in order; exit 1 if any failed."""
    env_ok = os.path.exists(os.path.join(HERE, "data", "isc_bulletin_events.csv"))
    fails = 0
    if env_ok:
        fails += bool(run(["build_combined_catalog.py"]))
    for cat in CATALOGS:
        for case in ([], ["--include-fixed"]):
            for s in SCRIPTS:
                fails += bool(run([s] + cat + case))
    fails += bool(run(["catalog_qc.py"]))
    for s in ("bvalue_robustness.py", "mc_through_time.py", "cluster_robustness.py"):
        fails += bool(run([s]))
    if "--baseline" in sys.argv:
        fails += bool(run(["reproduce_baseline.py", "--fig10"]))
    print("all done" if not fails else f"{fails} step(s) failed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
