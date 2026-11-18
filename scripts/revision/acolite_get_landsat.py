"""Download one Landsat C2 L1 product bundle via USGS M2M (credentials in ~/.usgs_m2m: 'username token')."""
import sys, os, time, tarfile, requests
pid, outdir = sys.argv[1], sys.argv[2]
U = "https://m2m.cr.usgs.gov/api/api/json/stable/"
user, tok = open(os.path.expanduser("~/.usgs_m2m")).read().split()
s = requests.Session()
def call(ep, body):
    r = s.post(U + ep, json=body, timeout=120); r.raise_for_status(); j = r.json()
    if j.get("errorCode"): raise SystemExit(f"{ep}: {j['errorCode']} {j.get('errorMessage')}")
    return j["data"]
key = call("login-token", {"username": user, "token": tok}); s.headers["X-Auth-Token"] = key
ds = "landsat_ot_c2_l1"
lid = "aco_" + pid[-30:]
call("scene-list-add", {"listId": lid, "datasetName": ds, "idField": "displayId", "entityId": pid})
opts = call("download-options", {"datasetName": ds, "listId": lid})
call("scene-list-remove", {"listId": lid})
o = [x for x in opts if x.get("available") and "Bundle" in x.get("productName", "")]
if not o: raise SystemExit("no bundle option " + pid)
req = call("download-request", {"downloads": [{"entityId": o[0]["entityId"], "productId": o[0]["id"]}], "label": lid})
urls = [d["url"] for d in req.get("availableDownloads", [])]
for _ in range(30):
    if urls: break
    time.sleep(20)
    ret = call("download-retrieve", {"label": lid})
    urls = [d["url"] for d in ret.get("available", [])]
if not urls: raise SystemExit("download not ready " + pid)
tar = os.path.join(outdir, pid + ".tar")
with s.get(urls[0], stream=True, timeout=600) as r:
    r.raise_for_status()
    with open(tar, "wb") as f:
        for c in r.iter_content(1 << 20): f.write(c)
d = os.path.join(outdir, pid); os.makedirs(d, exist_ok=True)
with tarfile.open(tar) as t:
    t.extractall(d, members=[m for m in t.getmembers() if not m.name.endswith(("B8.TIF",))], filter="data")
os.remove(tar)
call("logout", {})
print("OK", pid)
