"""Download one Sentinel-2 L1C SAFE from the Google public bucket (skips previews)."""
import sys, os, requests
pid, outdir = sys.argv[1], sys.argv[2]
tile = pid.split('_')[5][1:]
pre = f"tiles/{tile[:2]}/{tile[2]}/{tile[3:]}/{pid}.SAFE/"
B = "https://storage.googleapis.com/storage/v1/b/gcp-public-data-sentinel-2/o"
names, tok = [], None
while True:
    r = requests.get(B, params={"prefix": pre, "pageToken": tok} if tok else {"prefix": pre}, timeout=60).json()
    names += [o["name"] for o in r.get("items", [])]
    tok = r.get("nextPageToken")
    if not tok: break
if not names: sys.exit("NOT FOUND " + pid)
skip = ("_PVI.jp2", "_TCI.jp2", "$folder$")
for n in names:
    if n.endswith(skip) or n.endswith("/"): continue
    dst = os.path.join(outdir, n[len(pre) - len(pid) - 6:])
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with requests.get("https://storage.googleapis.com/gcp-public-data-sentinel-2/" + n, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(dst, "wb") as f:
            for c in r.iter_content(1 << 20): f.write(c)
for sub in ("AUX_DATA", "HTML", "rep_info"):
    os.makedirs(os.path.join(outdir, pid + ".SAFE", sub), exist_ok=True)
print("OK", pid, len(names))
