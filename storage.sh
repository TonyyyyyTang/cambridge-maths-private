#!/usr/bin/env bash
# Storage helper for the iCloud-backed raw binary store.
#
#   ./storage.sh link     (re)create private/raw and private/my-handwritten symlinks (new machine setup)
#   ./storage.sh check    list files iCloud has evicted to cloud-only placeholders ("dataless")
#   ./storage.sh fetch    force-download every evicted file, then re-check
#   ./storage.sh verify   compare sha256 of every file in manifest.csv; report missing/mismatched
#   ./storage.sh untracked  list files in raw/ that are not in manifest.csv
#
# Run `fetch` before any batch job that reads raw files, and `verify` after any batch download.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
STORE="$HOME/Library/Mobile Documents/com~apple~CloudDocs/Cambridge-Maths-Raw"
MANIFEST="$HERE/manifest.csv"
DIRS=(raw my-handwritten)

evicted() {
  for d in "${DIRS[@]}"; do
    find "$STORE/$d" -type f -flags +dataless 2>/dev/null || true
  done
}

case "${1:-}" in
  link)
    mkdir -p "$STORE/raw" "$STORE/my-handwritten"
    for d in "${DIRS[@]}"; do
      if [ -L "$HERE/$d" ]; then echo "ok: $d already linked"; continue; fi
      if [ -e "$HERE/$d" ]; then echo "ERROR: $HERE/$d exists and is not a symlink" >&2; exit 1; fi
      ln -s "$STORE/$d" "$HERE/$d" && echo "linked: $d -> $STORE/$d"
    done
    ;;
  check)
    n=$(evicted | wc -l | tr -d ' ')
    if [ "$n" = 0 ]; then echo "ok: no evicted files"; else evicted; echo "$n evicted file(s) — run ./storage.sh fetch"; exit 2; fi
    ;;
  fetch)
    evicted | while IFS= read -r f; do brctl download "$f"; done
    echo "requested downloads; waiting for iCloud..."
    for _ in $(seq 1 60); do
      [ "$(evicted | wc -l | tr -d ' ')" = 0 ] && { echo "ok: all files local"; exit 0; }
      sleep 5
    done
    echo "still evicted after 5 min:"; evicted; exit 2
    ;;
  verify)
    [ -f "$MANIFEST" ] || { echo "no manifest.csv yet"; exit 0; }
    bad=0
    # manifest columns: url,local_path,sha256,fetched_at,source  (local_path relative to private/)
    tail -n +2 "$MANIFEST" | while IFS=, read -r url path sha fetched source; do
      f="$HERE/$path"
      if [ ! -f "$f" ]; then echo "MISSING  $path"; bad=1; continue; fi
      if [ -n "$(find "$f" -flags +dataless)" ]; then echo "EVICTED  $path (run fetch first)"; bad=1; continue; fi
      actual=$(shasum -a 256 "$f" | cut -d' ' -f1)
      [ "$actual" = "$sha" ] || { echo "MISMATCH $path"; bad=1; }
    done
    echo "verify done"
    ;;
  untracked)
    [ -f "$MANIFEST" ] || { echo "no manifest.csv yet"; exit 0; }
    for d in "${DIRS[@]}"; do
      (cd "$HERE" && find "$d/" -type f ! -name '.DS_Store')
    done | while IFS= read -r p; do
      grep -qF ",$p," "$MANIFEST" || echo "UNTRACKED $p"
    done
    ;;
  *)
    sed -n '2,10p' "$0"; exit 1 ;;
esac
