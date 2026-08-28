"""Print the fitted membrane geometry so it can be compared across architectures."""
import platform, numpy, scipy
from memorient.barrel import fit_membrane
from memorient.contexts import get_context
from memorient.membrane import project_membrane, INTERFACE_WIDTH
from memorient.sasa import compute_sasa
from synthetic import make_barrel

print("machine:", platform.machine(), "| python", platform.python_version(),
      "| numpy", numpy.__version__, "| scipy", scipy.__version__)
try:
    cfg = numpy.__config__.show(mode="dicts")
    blas = (cfg or {}).get("Build Dependencies", {}).get("blas", {})
    print("blas:", blas.get("name"), blas.get("version"), "| openblas threads:",
          __import__("os").environ.get("OPENBLAS_NUM_THREADS", "unset"))
except Exception as exc:
    print("blas: could not introspect:", exc)
GN = get_context("gram_negative_om")
s = make_barrel(n_strands=12, strand_len=10, seed=0)
fit = fit_membrane(s, GN)
rsa = compute_sasa(s, n_points=120)["rsa"]
proj = project_membrane(s, fit, GN, ec_sign=1, rsa=rsa)
d = fit.half_thickness
depth = (s.ca - fit.centroid) @ fit.normal - fit.center
lim = d + INTERFACE_WIDTH
print(f"half_thickness d = {d:.6f}   limit = {lim:.6f}")
print(f"normal           = [{fit.normal[0]:+.6f} {fit.normal[1]:+.6f} {fit.normal[2]:+.6f}]")
print(f"depth min {depth.min():+.4f} (margin {(-depth.min())-lim:+.4f})  "
      f"max {depth.max():+.4f} (margin {depth.max()-lim:+.4f})")
print("zones:", sorted(set(proj.zone.tolist())))
