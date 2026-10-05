"""Embed saved nominal CSV rows and existing figures in an offline HTML viewer.

No ODE integration, fitting, interpolation or external network request occurs.
Run from any directory: python scripts/build_offline_viewer.py
"""
from pathlib import Path
import base64
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def build():
    datasets = []
    study = ROOT / "results/study_full_v1"
    for path in sorted(study.glob("*_nominal_trajectories.csv")):
        label = path.name.removesuffix("_nominal_trajectories.csv")
        model = label.split("_")[0]
        groups = {}
        with path.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                key = row["arm"] + " / " + row["policy"]
                groups.setdefault(key, []).append([float(row[k]) for k in ("t", "x", "s", "d")])
        datasets.append(dict(label=label, model=model, source=path.relative_to(ROOT).as_posix(),
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest(), series=groups))
    allowlist = json.loads((ROOT / "docs/source_allowlist.json").read_text(encoding="utf-8"))
    sources = [name for name in allowlist["copied_allowlist"]
               if name.startswith("results/observation_full_") and name.endswith("/latent_trajectories.csv")]
    if len(sources) != 1:
        raise ValueError("Exactly one frozen Stage D latent source is required")
    path = ROOT / sources[0]
    with path.open(encoding="utf-8", newline="") as stream:
        rows = [[float(row["time"]), float(row["x"]), None, float(row["d"])]
                for row in csv.DictReader(stream) if row["case"] == "nominal"]
    datasets.append(dict(label="M4_frozen_d015", model="M4", source=path.relative_to(ROOT).as_posix(),
                         sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         series={"both / maintenance (Stage D frozen d=0.15)": rows}))
    payload = json.dumps(datasets, separators=(",", ":"), allow_nan=False)
    figures = []
    captions = ("Figure 1. Saved nominal policies: maintenance and full withdrawal are distinct.",
                "Figure 2. Saved sampled-design target fractions; these are not patient probabilities.",
                "Figure 3. Saved history and rate diagnostics; synthetic, dimensionless examples.")
    for i, path in enumerate(sorted((ROOT / "results/manuscript_figures_v1/revised_layout").glob("figure_*.png"))):
        figures.append('<figure><img alt="' + captions[i] + '" src="data:image/png;base64,'
                       + base64.b64encode(path.read_bytes()).decode("ascii") + '"><figcaption>'
                       + captions[i] + '</figcaption></figure>')
    template = (ROOT / "docs/viewer_template.html").read_text(encoding="utf-8")
    output = template.replace("__SAVED_DATA__", payload).replace("__SAVED_FIGURES__", "\n".join(figures))
    (ROOT / "docs/model_viewer.html").write_text(output, encoding="utf-8")
    return {"datasets": len(datasets), "series": sum(len(d["series"]) for d in datasets),
            "rows": sum(len(s) for d in datasets for s in d["series"].values()),
            "figures": len(figures)}


if __name__ == "__main__":
    print(json.dumps(build(), sort_keys=True))
