"""X-066 support: weather flags for 2025, where the field station record (2020–2024) ends.

The laser runs to 2025-12 and the 2025 Sentinel-1 pairs exist (D-023), but the laser snow mask and the
wet / frozen overpass flags (X-050) need hourly air temperature, humidity and rain at the site. Source:
Open-Meteo's historical archive (ERA5-based reanalysis, free, no key) at the site centroid. Before it is
used, it is scored against the station on 2020–2024 with the station's own rules (``field.snow_mask``,
``field.surface_wetness_at``), and calibrated there: an air-temperature offset chosen to reproduce the station's
laser snow mask (the reanalysis is warmer than the peatland, most on cold nights), the mean air-temperature bias
for 'frozen at the overpass', and the humidity threshold that best reproduces the
station's wet flags. The same settings are applied to 2025. Outputs only to the hub (the laser is field data).

    PYTHONPATH=src python scripts/fusion_v2_flags_2025.py
"""
from __future__ import annotations

import io
import json
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from insar_wetlands import field
from insar_wetlands.config import load_config

HUB = field.field_root().parents[1]
RAW = HUB / "06_data" / "local_small" / "open_meteo"
OUT = HUB / "08_deliverables" / "fusion_v2" / "flags_2025"
HOUR = {"ascending": "16:36", "descending": "05:09"}
API = "https://archive-api.open-meteo.com/v1/archive"


def fetch(lat: float, lon: float, start: str, end: str) -> pd.DataFrame:
    q = urllib.parse.urlencode({"latitude": lat, "longitude": lon, "start_date": start, "end_date": end,
                                "hourly": "temperature_2m,relative_humidity_2m,precipitation",
                                "timezone": "GMT", "format": "csv"})
    with urllib.request.urlopen(f"{API}?{q}", timeout=300) as r:
        txt = r.read().decode()
    body = txt[txt.index("time,"):]                       # the CSV after the location header block
    d = pd.read_csv(io.StringIO(body))
    d["time"] = pd.to_datetime(d["time"], utc=True)
    return d.set_index("time").rename(columns={"temperature_2m (°C)": "Air_2m", "relative_humidity_2m (%)": "RH_2m",
                                               "precipitation (mm)": "Rain_mm_Tot"})


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    lon, lat = load_config()["site"]["centroid"]
    f = RAW / "open_meteo_rzecin_hourly_2020_2025.csv"
    if not f.exists():
        fetch(lat, lon, "2020-01-01", "2025-12-31").to_csv(f)
    om = pd.read_csv(f, parse_dates=["time"]).set_index("time")
    om.index = pd.to_datetime(om.index, utc=True)
    st = field.load_wtd_hourly()                               # station: Air_2m, RH_2m, Rain_mm_Tot (to 2024)

    # 1. agreement with the station on 2020–2024 (station stamps are at half past: interpolate onto them)
    omi = om.reindex(om.index.union(st.index)).interpolate("time").reindex(st.index)
    both = st[["Air_2m", "RH_2m"]].join(omi[["Air_2m", "RH_2m"]], rsuffix="_om", how="inner").dropna()
    agree = {"hours": len(both),
             "air_r": float(np.corrcoef(both.Air_2m, both.Air_2m_om)[0, 1]),
             "air_bias_c": float((both.Air_2m_om - both.Air_2m).mean()),
             "rh_r": float(np.corrcoef(both.RH_2m, both.RH_2m_om)[0, 1]),
             "rh_bias_pct": float((both.RH_2m_om - both.RH_2m).mean())}

    # 2. calibrate Open-Meteo to reproduce the station's flags on 2020–2024
    lz = field.load_laser()
    idx = lz.index[(lz.index >= "2021-06-01") & (lz.index < "2025-01-01")]
    s_st = field.snow_mask(st["Air_2m"], idx)
    best_dt = max(np.arange(0.0, 4.01, 0.25),
                  key=lambda dt: float((field.snow_mask(om["Air_2m"] - dt, idx) == s_st).mean()))
    om_adj = om.assign(Air_2m=om["Air_2m"] - best_dt)          # for the 72-h snow mask (coldest hours matter)
    om_ovp = om.assign(Air_2m=om["Air_2m"] - agree["air_bias_c"])  # for 'frozen at the overpass': the mean bias
    times = {t: [pd.Timestamp(f"{d.date()} {HOUR[t]}", tz="UTC") for d in pd.date_range("2020-01-01", "2024-12-31", freq="6D")]
             for t in HOUR}
    st_flags = {t: field.surface_wetness_at(st, tt) for t, tt in times.items()}
    def wet_agreement(rh):
        return np.mean([(field.surface_wetness_at(om_ovp, tt, rh_wet=rh).wet == st_flags[t].wet).mean()
                        for t, tt in times.items()])
    best_rh = max(np.arange(85.0, 99.1, 1.0), key=wet_agreement)
    agree.update({"air_offset_applied_c": float(best_dt), "rh_wet_threshold_open_meteo": float(best_rh),
                  "laser_snow72_agreement_2021_2024": float((field.snow_mask(om_adj["Air_2m"], idx) == s_st).mean()),
                  "laser_snow72_share_station": float(s_st.mean()),
                  "laser_snow72_share_open_meteo_calibrated": float(field.snow_mask(om_adj["Air_2m"], idx).mean())})
    rows = []
    for t, tt in times.items():
        a, b = st_flags[t], field.surface_wetness_at(om_ovp, tt, rh_wet=best_rh)
        for flag in ("wet", "frozen"):
            rows.append({"track": t, "flag": flag, "n": len(a), "station_share": float(a[flag].mean()),
                         "open_meteo_share": float(b[flag].mean()), "agreement": float((a[flag] == b[flag]).mean())})
    flags_cmp = pd.DataFrame(rows)

    # 3. 2025: laser snow mask and overpass flags from the calibrated Open-Meteo
    i25 = lz.index[lz.index.year == 2025]
    snow25 = field.snow_mask(om_adj["Air_2m"], i25)
    pd.DataFrame({"time_utc": i25, "snow_72h_open_meteo": snow25.values}).to_csv(OUT / "laser_snow72_2025.csv", index=False)
    ov = []
    for t in HOUR:
        dates = sorted({pd.Timestamp(x) for p in (HUB / "05_code" / "local" / "s1_2025" / f"hyp3_cropped_{t}").iterdir()
                        if p.is_dir() for x in p.name.split("_")})
        w = field.surface_wetness_at(om_ovp, [pd.Timestamp(f"{d.date()} {HOUR[t]}", tz="UTC") for d in dates], rh_wet=best_rh)
        ov.append(w.assign(date=[d.date() for d in dates], track=t))
    pd.concat(ov).to_csv(OUT / "wetness_at_overpasses_2025_open_meteo.csv", index=False)

    flags_cmp.to_csv(OUT / "open_meteo_vs_station_flags.csv", index=False)
    (OUT / "open_meteo_vs_station.json").write_text(json.dumps(agree, indent=2))
    print(json.dumps(agree, indent=2))
    print(flags_cmp.round(3).to_string(index=False))
    print("2025 laser hours usable for snow:", int((~snow25).sum()), "of", len(i25))


if __name__ == "__main__":
    main()
