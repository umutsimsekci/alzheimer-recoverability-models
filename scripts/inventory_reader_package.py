"""Refresh public reader hashes after intentional documentation/metadata edits.

This preserves frozen source provenance and refuses changes to exact science
copies. It does not add Git/runtime/environment files to the public inventory.
"""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_NAMES = {"MANIFEST.json", "SHA256SUMS", "FILE_INVENTORY.csv"}


def refresh():
    provenance = json.loads((ROOT / "docs/source_provenance.json").read_text(encoding="utf-8"))
    for source in provenance["files"]:
        path = ROOT / source["path"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != source["reader_copy_sha256"]:
            raise ValueError("Frozen reader source changed: " + source["path"])
    inference_provenance = ROOT / "docs/inference_source_provenance.json"
    if inference_provenance.exists():
        for source in json.loads(inference_provenance.read_text(encoding="utf-8"))["files"]:
            path = ROOT / source["path"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != source["reader_copy_sha256"]:
                raise ValueError("Frozen inference source changed: " + source["path"])
    paths = []
    for path in sorted(ROOT.rglob("*")):
        relative = path.relative_to(ROOT)
        if not path.is_file() or path.is_symlink():
            continue
        if any(part in {".git", ".venv", "__pycache__", ".pytest_cache"} for part in relative.parts):
            continue
        if relative.parts[0] == "results" and len(relative.parts) > 1 and relative.parts[1].startswith("reproduced"):
            continue
        if path.suffix in {".pyc", ".pyo"} or path.name.endswith(".local.json") or path.name == ".env":
            continue
        if len(relative.parts) == 1 and path.name in EXCLUDED_NAMES:
            continue
        paths.append(path)
    entries = [dict(path=path.relative_to(ROOT).as_posix(), bytes=path.stat().st_size,
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in paths]
    with (ROOT / "FILE_INVENTORY.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(entries)
    path = ROOT / "FILE_INVENTORY.csv"
    entries.append(dict(path=path.name, bytes=path.stat().st_size,
                        sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    manifest = dict(schema_version=1, evidence_class="synthetic_simulation_and_technical_verification",
                    release_status="reader code and synthetic data package v2.0.0; no archival DOI or reuse license assigned",
                    public_repository_url="https://github.com/umutsimsekci/alzheimer-recoverability-models",
                    source_snapshot="original dynamics plus revision_20261019 inference allowlist",
                    provenance_file="docs/source_provenance.json", inference_provenance_file="docs/inference_source_provenance.json", files=entries,
                    self_hash_policy="Manifest excludes itself and SHA256SUMS. CSV excludes itself and both hash-control files.")
    (ROOT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    checksums = [(e["path"], e["sha256"]) for e in entries]
    checksums.append(("MANIFEST.json", hashlib.sha256((ROOT / "MANIFEST.json").read_bytes()).hexdigest()))
    (ROOT / "SHA256SUMS").write_text("".join(digest+"  "+name+"\n" for name, digest in sorted(checksums)), encoding="utf-8")
    return {"content_files": len(paths), "manifest_entries": len(entries), "checksum_entries": len(checksums)}


if __name__ == "__main__":
    print(json.dumps(refresh(), sort_keys=True))
