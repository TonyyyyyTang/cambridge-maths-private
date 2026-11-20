# cambridge-maths-private

Private layer for Cambridge All-in-one Maths. **Never publish.**

- Text (transcriptions, links.yaml, marking notes, manifest.csv) is tracked here.
- Binaries (PDFs, scans) live in iCloud Drive at `Cambridge-Maths-Raw/`, symlinked as `raw/` and `my-handwritten/`.
- `manifest.csv`: one row per binary — `url,local_path,sha256,fetched_at,source` (`local_path` relative to this folder, e.g. `raw/tripos/IB/2019/P1.pdf`).

## Storage rules
1. New machine: `./storage.sh link`.
2. Before a batch job reads raw files: `./storage.sh fetch` (iCloud may evict files to cloud-only placeholders).
3. After a batch download: add rows to `manifest.csv`, then `./storage.sh verify` and `./storage.sh untracked`.
4. Never create a git repo inside iCloud Drive. Never move these repos into iCloud.
5. Filenames: ASCII only, `[A-Za-z0-9._-]`, no spaces or colons. iCloud conflict copies ("file 2.pdf") are duplicates — delete after checking sha256.
