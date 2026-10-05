"""Read-only integrity and saved-viewer provenance checks; no simulations."""
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import sys

ROOT = Path(__file__).resolve().parents[1]


def verify():
    manifest = json.loads((ROOT / "MANIFEST.json").read_text(encoding="utf-8"))
    failures = []
    for entry in manifest["files"]:
        path = ROOT / entry["path"]
        if not path.is_file() or path.is_symlink():
            failures.append("Missing or symlink: " + entry["path"])
            continue
        if path.stat().st_size != entry["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            failures.append("Hash/size mismatch: " + entry["path"])
    checksum_entries = 0
    for line in (ROOT / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        path = ROOT / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            failures.append("SHA256SUMS mismatch: " + name)
        checksum_entries += 1
    html = (ROOT / "docs/model_viewer.html").read_text(encoding="utf-8")
    match = re.search(r'<script id="saved-data" type="application/json">(.*?)</script>', html, re.S)
    data = json.loads(match.group(1))
    count = 0
    for dataset in data:
        path = ROOT / dataset["source"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != dataset["sha256"]:
            failures.append("Viewer source hash mismatch: " + dataset["source"])
        expected = {}
        with path.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                if dataset["model"] == "M4":
                    if row["case"] != "nominal":
                        continue
                    key = "both / maintenance (Stage D frozen d=0.15)"
                    values = [float(row["time"]), float(row["x"]), None, float(row["d"])]
                else:
                    key = row["arm"] + " / " + row["policy"]
                    values = [float(row[k]) for k in ("t", "x", "s", "d")]
                expected.setdefault(key, []).append(values)
        if expected != dataset["series"]:
            failures.append("Embedded rows differ from CSV: " + dataset["label"])
        for values in dataset["series"].values():
            if any(not math.isfinite(v) for row in values for v in row if v is not None):
                failures.append("Nonfinite viewer row: " + dataset["label"])
            if any(b[0] <= a[0] for a, b in zip(values, values[1:])):
                failures.append("Non-increasing saved times: " + dataset["label"])
            count += len(values)
    if len(data) != 7 or sum(len(d["series"]) for d in data) != 37:
        failures.append("Unexpected viewer design coverage")
    if len(re.findall(r'<img[^>]*src="data:image/png;base64,', html)) != 3:
        failures.append("Expected three embedded existing figures")
    if re.search(r'<(?:script|link|iframe)[^>]*(?:src|href)="https?://|\bfetch\s*\(|XMLHttpRequest', html):
        failures.append("External network dependency found in viewer")
    return {"manifest_entries": len(manifest["files"]), "checksum_entries": checksum_entries,
            "viewer_datasets": len(data), "viewer_series": sum(len(d["series"]) for d in data),
            "viewer_exact_saved_rows": count, "embedded_existing_figures": 3,
            "failures": failures}


if __name__ == "__main__":
    result = verify()
    print(json.dumps(result, indent=2))
    sys.exit(bool(result["failures"]))
