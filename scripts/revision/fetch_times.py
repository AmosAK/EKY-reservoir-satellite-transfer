"""Download WQP activity start times for the training stations (POST, 40 sites per request)."""
import io, sys, time
import pandas as pd, requests

st = pd.read_csv(sys.argv[1]).station.tolist()
out_path = sys.argv[2]
URL = "https://www.waterqualitydata.us/data/Activity/search?mimeType=csv&zip=no"
cols = ["MonitoringLocationIdentifier", "ActivityStartDate", "ActivityStartTime/Time",
        "ActivityStartTime/TimeZoneCode", "ActivityTypeCode", "ActivityMediaName"]
parts, failed = [], []
chunks = [st[i:i + 40] for i in range(0, len(st), 40)]
for i, ch in enumerate(chunks):
    body = {"siteid": ch, "startDateLo": "01-01-2013", "startDateHi": "12-31-2026"}
    for attempt in range(3):
        try:
            r = requests.post(URL, json=body, timeout=300)
            r.raise_for_status()
            d = pd.read_csv(io.StringIO(r.text), usecols=lambda c: c in cols, dtype=str)
            parts.append(d); break
        except Exception as e:
            if attempt == 2: failed.append(i); print("chunk", i, "failed:", e, flush=True)
            time.sleep(5)
    print(f"chunk {i + 1}/{len(chunks)} rows so far {sum(len(p) for p in parts)}", flush=True)
    time.sleep(0.5)
out = pd.concat(parts).drop_duplicates()
out.columns = [c.replace("ActivityStartTime/", "").replace("MonitoringLocationIdentifier", "station")
               .replace("ActivityStartDate", "date") for c in out.columns]
out.to_csv(out_path, index=False)
print("DONE rows", len(out), "stations", out.station.nunique(), "with time", out["Time"].notna().sum(), "failed chunks", failed)
