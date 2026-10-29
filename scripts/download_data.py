"""Download the Kang et al. 2018 PBMC dataset (h5ad) and verify its checksum.

The file is the same one served by pertpy's `pt.dt.kang_2018()` loader; fetching it
directly avoids installing pertpy just for a download.
"""
import argparse
import hashlib
import urllib.request
from pathlib import Path


def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--sha256", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part")
    print(f"Downloading {args.url}")
    urllib.request.urlretrieve(args.url, tmp)
    digest = sha256sum(tmp)
    if digest != args.sha256:
        tmp.unlink()
        raise SystemExit(f"Checksum mismatch: got {digest}, expected {args.sha256}")
    tmp.rename(out)
    print(f"Saved {out} ({out.stat().st_size / 1e6:.1f} MB), sha256 OK")


if __name__ == "__main__":
    main()
