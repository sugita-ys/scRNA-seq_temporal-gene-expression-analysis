from pathlib import Path
import requests
import sys

URL = "https://allen-developmental-mouse-atlas.s3.amazonaws.com/scRNA/DevVIS_scRNA_processed.h5ad"

def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/DevVIS_scRNA_processed.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(URL, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done = 0
        chunk = 1024 * 1024

        with open(out, "wb") as f:
            for part in r.iter_content(chunk_size=chunk):
                if part:
                    f.write(part)
                    done += len(part)
                    if total:
                        pct = 100 * done / total
                        print(f"\r{pct:6.2f}%  {done/(1024*1024):.1f}MB/{total/(1024*1024):.1f}MB", end="")
    print(f"\nSaved: {out}")

if __name__ == "__main__":
    main()
