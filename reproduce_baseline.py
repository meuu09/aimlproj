"""Run the original Him_ides scripts unchanged on the bundled catalog; outputs go to results_baseline/."""
import glob
import logging
import os
import runpy
import shutil
import subprocess
import sys
import types

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from seis_utils import max_curv_eqcorrscan

HERE = os.path.dirname(os.path.abspath(__file__))
UPSTREAM, COMMIT = "https://github.com/kemichai/Him_ides", "33acd38"
ORIG = os.path.join(HERE, "Him_ides", "python")
OUT = os.path.join(HERE, "results_baseline")
os.makedirs(OUT, exist_ok=True)


def fetch_upstream():
    """Clone the original repository at the commit we reproduced (it is not redistributed with this project)."""
    if os.path.isdir(ORIG):
        return
    dest = os.path.dirname(ORIG)
    subprocess.run(["git", "clone", "--quiet", UPSTREAM, dest], check=True)
    subprocess.run(["git", "-C", dest, "checkout", "--quiet", COMMIT], check=True)


fetch_upstream()
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


def calc_max_curv(magnitudes, bin_size=0.5, plotvar=False):
    """Stand-in for eqcorrscan 0.4.2 calc_max_curv (eqcorrscan no longer builds)."""
    return max_curv_eqcorrscan(magnitudes, bin_size)


def calc_b_value(*a, **k):
    """Imported by temporal_plots.py but never called."""
    raise NotImplementedError("not used by the Him_ides scripts")


pkg = types.ModuleType("eqcorrscan"); utils = types.ModuleType("eqcorrscan.utils")
mc = types.ModuleType("eqcorrscan.utils.mag_calc")
mc.calc_max_curv, mc.calc_b_value = calc_max_curv, calc_b_value
pkg.utils, utils.mag_calc = utils, mc
sys.modules.update({"eqcorrscan": pkg, "eqcorrscan.utils": utils, "eqcorrscan.utils.mag_calc": mc})

_shown = [0]


def _show(*a, **k):
    """Headless plt.show(): save the current figure instead of displaying it."""
    _shown[0] += 1
    plt.savefig(os.path.join(OUT, f"_shown_{_shown[0]}.png"), dpi=200, bbox_inches="tight")
    plt.close("all")


plt.show = _show

sys.path.insert(0, ORIG)
os.chdir(ORIG)
before = set(os.listdir(ORIG))

import temporal_plots as tp  # noqa: E402  (imports our eqcorrscan stand-in)
fits = []
_orig_calc_bvalues = tp.calc_bvalues


def _logged_calc_bvalues(mag, Mc, M_max, mag_steps, cdf):
    """Original calc_bvalues, also recording N, Mc and b for baseline_bvalues.csv."""
    fit_fn, fit = _orig_calc_bvalues(mag, Mc, M_max, mag_steps, cdf)
    used = [m for m in mag if m >= Mc + 0.4]
    fits.append(dict(N=len(mag), Mc_maxcurv=Mc, min_mag_in_fit=min(used), b=abs(fit[0]), mags=list(mag)))
    return fit_fn, fit


tp.calc_bvalues = _logged_calc_bvalues


def _axes3d_attached(fig, *args, **kwargs):
    """Axes3D(fig) no longer attaches the axes to the figure (matplotlib >= 3.4),
    which crashed Figure_10.py; create the 3-D axes the current way instead."""
    return fig.add_subplot(111, projection="3d")


try:
    print("Running Him_ides/python/Figures_5-9.py ...")
    runpy.run_path(os.path.join(ORIG, "Figures_5-9.py"), run_name="__main__")
    n_shown_59 = _shown[0]
    if "--fig10" in sys.argv:
        import mpl_toolkits.mplot3d as _m3d
        _m3d.Axes3D = _axes3d_attached
        print("Running Him_ides/python/Figure_10.py ...")
        runpy.run_path(os.path.join(ORIG, "Figure_10.py"), run_name="__main__")
finally:
    for f in sorted(set(os.listdir(ORIG)) - before):
        if f.lower().endswith((".png", ".jpg", ".pdf")):
            shutil.move(os.path.join(ORIG, f), os.path.join(OUT, f))

for f in glob.glob(os.path.join(OUT, "_shown_*.png")):
    k = int(os.path.basename(f)[7:-4])
    if k > n_shown_59:
        os.replace(f, os.path.join(OUT, "Fig_10.png" if k == n_shown_59 + 1 else f"Fig_10_{k - n_shown_59}.png"))
    else:
        os.remove(f)

if os.path.exists(os.path.join(OUT, "Fig_7.png")):
    os.replace(os.path.join(OUT, "Fig_7.png"), os.path.join(OUT, "Fig_7.jpg"))

names = ["Full catalog", "Cluster 1", "Cluster 2"]
paper = {"Full catalog": "Mc 2.4, b 0.91", "Cluster 1": "Mc 2.85, b 0.96", "Cluster 2": "not reported (61 events)"}
lines = ["group,N,Mc_maxcurv_raw,Mc_released_code,min_mag_in_fit,b_lsq,paper"]
mag_rows = ["group,mag"]
for n, f in zip(names, fits):
    lines.append(f"{n},{f['N']},{f['Mc_maxcurv']},{f['Mc_maxcurv'] + 0.4:.2f},{f['min_mag_in_fit']},{f['b']:.3f},"
                 f'"{paper[n]}"')
    mag_rows += [f"{n},{m}" for m in f["mags"]]
open(os.path.join(OUT, "baseline_magnitudes.csv"), "w").write("\n".join(mag_rows) + "\n")
open(os.path.join(OUT, "baseline_bvalues.csv"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("Outputs in", OUT, ":", sorted(os.listdir(OUT)))
