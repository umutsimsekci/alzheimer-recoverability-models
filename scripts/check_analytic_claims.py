"""Independently verify polynomial identities with SymPy, without importing model.py."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as sp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = (root / args.output).resolve()
    if not output.is_relative_to(root / "results") or output.exists():
        raise ValueError("Use a new JSON output inside results")
    started = time.perf_counter()
    x, a, b, r, k = sp.symbols("x a b r k", real=True)
    f = a * (1-x) + b*x*x*(1-x) - r*x
    substitutions = {a: b*x*x*(1-2*x), r: 2*b*x*(1-x)**2}
    nominal = f.subs({a: sp.Rational(6,100), b: 4, r: sp.Rational(105,100)})
    roots = sp.solve(nominal, x)
    checks = {
        "fold_drift_identity": sp.simplify(f.subs(substitutions)) == 0,
        "fold_derivative_identity": sp.simplify(sp.diff(f,x).subs(substitutions)) == 0,
        "cusp_second_derivative": sp.simplify(sp.diff(f,x,2).subs(x, sp.Rational(1,3))) == 0,
        "cusp_drive": sp.simplify(substitutions[a].subs(x, sp.Rational(1,3))-b/27) == 0,
        "cusp_repair": sp.simplify(substitutions[r].subs(x, sp.Rational(1,3))-8*b/27) == 0,
        "common_rate_scaling": sp.expand(f.subs({a:k*a,b:k*b,r:k*r}, simultaneous=True)-k*f) == 0,
        "nominal_separatrix": sp.simplify(nominal.subs(x,sp.Rational(2,5))) == 0,
        "repair_ordering_derivative": sp.diff(f,r) == -x,
        "drive_ordering_derivative": sp.diff(f,a) == 1-x,
    }
    root_rows = [{"x":float(v), "slope":float(sp.diff(nominal,x).subs(x,v))} for v in roots]
    paths = ["protocol/study_design_v1.md", "protocol/analytical_results.md", "scripts/check_analytic_claims.py"]
    result = {"evidence_class":"symbolic_mathematical_verification", "empirical_calibration":False,
              "checks":checks, "all_passed":all(checks.values()), "nominal_roots":sorted(root_rows,key=lambda row:row["x"]),
              "sympy_version":sp.__version__, "elapsed_seconds":time.perf_counter()-started,
              "input_sha256":{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in paths}}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"all_passed":result["all_passed"],"checks":len(checks),"output":str(output.relative_to(root))}))
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
