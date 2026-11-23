"""
Baseline reproduction: run the ORIGINAL Him_ides scripts, unchanged, on the
authors' bundled catalog and collect their figures in results_baseline/.

Why a wrapper instead of running `python Figures_5-9.py` directly:
  * the original code imports eqcorrscan (only for `calc_max_curv`), which no
    longer builds on current Python. We inject the function body copied verbatim
    from eqcorrscan 0.4.2 (the version the authors pinned, see
    Him_ides/python/README.md), so the numbers are computed with the same code.
  * the scripts call plt.show() and write PNGs into the current folder; we run
    them headless from Him_ides/python and move the outputs afterwards.
  * we also record the Mc / b-value numbers that the original code prints only
    inside figure legends, so they can be compared with the paper.

Run:  .venv/bin/python reproduce_baseline.py          (Figures 5-9)
      .venv/bin/python reproduce_baseline.py --fig10  (also Figure 10, slower)
"""
import glob
import logging
import os
import runpy
import shutil
import sys
import types

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(HERE, "Him_ides", "python")
OUT = os.path.join(HERE, "results_baseline")
os.makedirs(OUT, exist_ok=True)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)  # 'family: normal' spam


# ---- eqcorrscan 0.4.2 stand-in (calc_max_curv copied verbatim, plotting removed) ----
def calc_max_curv(magnitudes, bin_size=0.5, plotvar=False):
    """eqcorrscan 0.4.2 eqcorrscan/utils/mag_calc.py::calc_max_curv (Chamberlain et al. 2018)."""
    min_bin, max_bin = int(min(magnitudes)), int(max(magnitudes) + 1)
    bins = np.arange(min_bin, max_bin + bin_size, bin_size)
    df, bins = np.histogram(magnitudes, bins)
    grad = (df[1:] - df[0:-1]) / bin_size
    curvature = (grad[1:] - grad[0:-1]) / bin_size
    max_curv = bins[np.argmax(np.abs(curvature))] + bin_size
    return float(max_curv)


def calc_b_value(*a, **k):  # imported by temporal_plots.py but never called
    raise NotImplementedError("not used by the Him_ides scripts")


pkg = types.ModuleType("eqcorrscan"); utils = types.ModuleType("eqcorrscan.utils")
mc = types.ModuleType("eqcorrscan.utils.mag_calc")
mc.calc_max_curv, mc.calc_b_value = calc_max_curv, calc_b_value
pkg.utils, utils.mag_calc = utils, mc
sys.modules.update({"eqcorrscan": pkg, "eqcorrscan.utils": utils, "eqcorrscan.utils.mag_calc": mc})

# ---- headless: plt.show() saves any figure that was not saved explicitly ----
_shown = [0]
def _show(*a, **k):
    _shown[0] += 1
    plt.savefig(os.path.join(OUT, f"_shown_{_shown[0]}.png"), dpi=200, bbox_inches="tight")
    plt.close("all")
plt.show = _show

# ---- run the original scripts from their own folder ----
sys.path.insert(0, ORIG)
os.chdir(ORIG)
before = set(os.listdir(ORIG))

import temporal_plots as tp  # noqa: E402  (imports our eqcorrscan stand-in)
fits = []
_orig_calc_bvalues = tp.calc_bvalues
def _logged_calc_bvalues(mag, Mc, M_max, mag_steps, cdf):
    fit_fn, fit = _orig_calc_bvalues(mag, Mc, M_max, mag_steps, cdf)
    used = [m for m in mag if m >= Mc + 0.4]
    fits.append(dict(N=len(mag), Mc_maxcurv=Mc, Mc_label=min(used), b=abs(fit[0])))
    return fit_fn, fit
tp.calc_bvalues = _logged_calc_bvalues

print("Running Him_ides/python/Figures_5-9.py ...")
runpy.run_path(os.path.join(ORIG, "Figures_5-9.py"), run_name="__main__")
if "--fig10" in sys.argv:
    print("Running Him_ides/python/Figure_10.py ...")
    runpy.run_path(os.path.join(ORIG, "Figure_10.py"), run_name="__main__")

# ---- move new outputs out of the vendored folder ----
for f in sorted(set(os.listdir(ORIG)) - before):
    if f.lower().endswith((".png", ".jpg", ".pdf")):
        shutil.move(os.path.join(ORIG, f), os.path.join(OUT, f))
for f in glob.glob(os.path.join(OUT, "_shown_*.png")):  # figures already saved by the scripts
    if "--fig10" not in sys.argv:
        os.remove(f)

# Fig_7 is written as JPEG with a .png name by the original code; keep a correct extension
if os.path.exists(os.path.join(OUT, "Fig_7.png")):
    os.replace(os.path.join(OUT, "Fig_7.png"), os.path.join(OUT, "Fig_7.jpg"))

names = ["Full catalog", "Cluster 1", "Cluster 2"]  # order of calls inside multi_mfd_plot
paper = {"Full catalog": "Mc 2.4, b 0.91", "Cluster 1": "Mc 2.85, b 0.96", "Cluster 2": "not reported (61 events)"}
lines = ["group,N,Mc_maxcurv_raw,Mc_plus0.4_label,b_lsq,paper"]
for n, f in zip(names, fits):
    lines.append(f"{n},{f['N']},{f['Mc_maxcurv']},{f['Mc_label']},{f['b']:.3f},{paper[n]}")
open(os.path.join(OUT, "baseline_bvalues.csv"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("Outputs in", OUT, ":", sorted(os.listdir(OUT)))
