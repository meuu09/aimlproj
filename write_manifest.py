"""Write results/MANIFEST.json: input checksums, software versions, seeds, git commit and output checksums."""
import hashlib
import json
import os
import platform
import subprocess
import time
from importlib import metadata

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = {"seis_utils.b_bootstrap": 0, "seis_utils.mc_sliding": 0, "seis_utils.cv_poisson_test": 0,
         "seis_utils.migration_fit permutation": 0, "cluster_robustness Monte Carlo": 0}
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "scikit-learn", "shapely", "obspy", "pytest"]


def sha(path):
    """SHA-256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tree(root, exts):
    """Checksums of every file under root with one of the extensions, keyed by relative path."""
    out = {}
    for d, _, files in os.walk(os.path.join(HERE, root)):
        for f in sorted(files):
            if f.endswith(exts):
                p = os.path.join(d, f)
                out[os.path.relpath(p, HERE)] = sha(p)
    return dict(sorted(out.items()))


def main():
    """Collect and write the manifest."""
    git = lambda *a: subprocess.run(["git", *a], cwd=HERE, capture_output=True, text=True).stdout.strip()
    versions = {}
    for p in PACKAGES:
        try:
            versions[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            versions[p] = None
    m = dict(date=time.strftime("%Y-%m-%dT%H:%M:%S%z"), git_commit=git("rev-parse", "HEAD"),
             git_dirty=bool(git("status", "--porcelain", "--", "*.py")), python=platform.python_version(),
             platform=platform.platform(), packages=versions, seeds=SEEDS,
             inputs=tree("data", (".csv", ".txt", ".geojson")),
             outputs=tree("results", (".csv", ".md")) | tree("results_baseline", (".csv", ".md")))
    json.dump(m, open(os.path.join(HERE, "results", "MANIFEST.json"), "w"), indent=1)
    print(f"{len(m['inputs'])} inputs, {len(m['outputs'])} outputs -> results/MANIFEST.json")


if __name__ == "__main__":
    main()
