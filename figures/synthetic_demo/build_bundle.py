#!/usr/bin/env python3
"""Package current, verified synthetic figures with their source and audit evidence."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
OUTPUT = ROOT.parent / "lipohypolep_synthetic_figure_pack.zip"


def main():
    subprocess.run([sys.executable, str(ROOT / "validate_figures.py"), "--check-report"], check=True)
    # Honor .gitignore and restrict the bundle to project source/demo material.
    # Experimental CSVs are never added by this packaging helper.
    roots = ("README.md", ".gitignore", "figures", "analysis", "tests", "docs")
    try:
        checkout = subprocess.check_output(["git", "rev-parse", "--show-toplevel"],
                                           cwd=REPO, stderr=subprocess.DEVNULL)
        if Path(checkout.decode().strip()).resolve() != REPO:
            raise ValueError("The source directory is not a Git checkout")
        listing = subprocess.check_output([
            "git", "ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", *roots,
        ], cwd=REPO)
        names = sorted({name.decode("utf-8") for name in listing.split(b"\0") if name})
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError):
        # The downloadable figure pack intentionally has no .git directory.
        names = sorted({str(path.relative_to(REPO)) for name in roots
                        for path in ([REPO / name] if (REPO / name).is_file()
                                     else (REPO / name).rglob("*")) if path.is_file()})
    entries = {}
    for name in names:
        path = REPO / name
        if name != ".gitignore" and any(part.startswith(".") for part in Path(name).parts):
            continue
        if path == OUTPUT or path.suffix in {".zip", ".csv", ".pyc"} or "__pycache__" in path.parts:
            continue
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Bundle source must be an ordinary file: {name}")
        if name.startswith("analysis/") and path.suffix not in {".py", ".md", ".txt"}:
            continue
        entries[name] = path.read_bytes()
    manifest = {"scope": "Synthetic figure development only; no experimental CSVs.",
                "files": {name: {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
                          for name, content in entries.items()}}
    entries["BUNDLE_MANIFEST.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, content in sorted(entries.items()):
            member = zipfile.ZipInfo(name, date_time=(2000, 1, 1, 0, 0, 0))
            member.compress_type = zipfile.ZIP_DEFLATED
            member.external_attr = 0o644 << 16
            archive.writestr(member, content)
    payload = buffer.getvalue()
    if len(payload) >= 40_000_000:
        raise ValueError("Figure pack exceeds the 40 MB delivery limit")
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        if archive.testzip() is not None:
            raise ValueError("Figure pack failed CRC validation")
        for name, details in manifest["files"].items():
            if hashlib.sha256(archive.read(name)).hexdigest() != details["sha256"]:
                raise ValueError(f"Figure pack content mismatch: {name}")
    with tempfile.NamedTemporaryFile(dir=OUTPUT.parent, suffix=".zip", delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(payload)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        temporary.replace(OUTPUT)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Saved {OUTPUT.name}: {len(entries)} files, {len(payload):,} bytes; CRC and SHA-256 verified")


if __name__ == "__main__":
    main()
