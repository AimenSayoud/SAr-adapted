"""Candidate Sentinel-1 burst pairs for extending the record back to 2020–2021 (D-020 / X-050).

Read-only: queries ASF's public search API (no credentials, no credits) for the two bursts in
config.yaml, and builds the small-baseline network with the pipeline's own
``hyp3.jobs.build_sbas_pairs`` (S1A + S1B: 6-day repeat until S1B failed in Dec 2021), plus the
bridge pairs into the first 2022 dates so the new network joins the existing one. Submitting the
pairs is a separate, confirmed step (``hyp3.jobs.submit_pairs`` needs its confirmation phrase).

    PYTHONPATH=src python scripts/s1_extension_candidates.py --out <dir>
"""
from __future__ import annotations

import argparse
import io
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

from insar_wetlands.hyp3.jobs import build_sbas_pairs

REPO = Path(__file__).resolve().parents[1]
ASF = "https://api.daac.asf.alaska.edu/services/search/param"


def search_bursts(full_burst_id: str, start: str, end: str) -> pd.DataFrame:
    q = urllib.parse.urlencode({"dataset": "SLC-BURST", "fullBurstID": full_burst_id, "polarization": "VV",
                                "start": start, "end": end, "output": "csv"})
    with urllib.request.urlopen(f"{ASF}?{q}", timeout=180) as r:
        d = pd.read_csv(io.BytesIO(r.read()))
    d["date"] = pd.to_datetime(d["Start Time"]).dt.tz_localize(None).dt.normalize()
    return d.drop_duplicates("date").sort_values("date")[["Granule Name", "Platform", "Start Time", "date"]]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", default="2020-01-01")
    ap.add_argument("--end", default="2022-01-01")
    ap.add_argument("--max-dt", type=int, default=24)
    ap.add_argument("--max-forward", type=int, default=4)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tracks = yaml.safe_load((REPO / "config" / "config.yaml").read_text())["sentinel1"]["tracks"]
    summary = []
    for track, t in tracks.items():
        d = search_bursts(t["burst_id"], f"{a.start}T00:00:00Z", "2022-03-01T00:00:00Z")
        pre = d[d.date < a.end]
        bridge = d[d.date >= a.end].head(4)
        pr = build_sbas_pairs(pd.concat([pre.date, bridge.date]), max_temporal_days=a.max_dt,
                              max_pairs_per_date=a.max_forward)
        pr = pr[pr.ref_date < a.end]
        pr.to_csv(out / f"candidate_pairs_{track}.csv", index=False)
        pre.to_csv(out / f"acquisitions_{track}.csv", index=False)
        summary.append({"track": track, "burst": t["burst_id"], "acquisitions": len(pre),
                        **{f"acq_{k}": int(v) for k, v in pre.Platform.value_counts().items()},
                        "pairs": len(pr), "bridge_pairs_into_2022": int((pr.sec_date >= a.end).sum()),
                        "max_dt_days": a.max_dt, "max_forward_per_date": a.max_forward})
    s = pd.DataFrame(summary)
    s.to_csv(out / "summary.csv", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
