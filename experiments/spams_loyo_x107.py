"""X-107 amendment_3 — the anchor's own spread: SPAMS refitted on the design years, leaving one year out each time.

The re-levelled estimators (v4, v4n) set each winter's level step from SPAMS, fitted once on the design years
(X-106). Their error statement (X-107 amendment_3) is the spread of that step over five refits, each without one
design year, with the freeze's settings. The five parameter sets are computed here, from design years only, and written
to ``08_deliverables/insar_full_test/v4/spams_loyo_params.yaml`` before any held-out laser year is delivered; X-107 then
evaluates them on the held-out winters. A check first reproduces the frozen fit (all design years) exactly.

Reads the as-delivered cube snapshot (as the freeze did: X-107 amendment) and the hub; writes only to the hub.

    RZECIN_CUBE_ROOT=../../06_data/cube_snapshots/as_delivered_2026-10-15 make run E=spams_loyo_x107
"""
from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import freeze_v4_x106 as F  # noqa: E402  (DESIGN, LAG, the frozen parameter file)
import spams_baseline_x112 as SB  # noqa: E402  (the daily weather/laser table the freeze fitted on)

from insar_wetlands import spams  # noqa: E402

OUT = F.V4DIR / "spams_loyo_params.yaml"


def main() -> None:
    frozen = yaml.safe_load(F.PARAMS.read_text(encoding="utf-8"))["spams"]
    d = SB.table()
    P, E, dh = d.rain_era5_mm.to_numpy(), d.et0_mm.to_numpy(), d.dh_mm.to_numpy()
    base = (d.dh_mm.notna() & ~d.frozen_win & d.year.isin(F.DESIGN)).to_numpy()
    full = spams.fit(P, E, dh, base, F.LAG)
    same = (np.isclose(full["x_p"], frozen["x_p"], rtol=0, atol=1e-12) and
            np.isclose(full["x_e"], frozen["x_e"], rtol=0, atol=1e-12) and int(full["tau"]) == int(frozen["tau_days"]))
    if not same:
        sys.exit(f"stop: the all-years fit does not reproduce the frozen SPAMS ({full} vs {frozen})")
    fits = {}
    for y in F.DESIGN:
        par = spams.fit(P, E, dh, base & (d.year != y).to_numpy(), F.LAG)
        fits[str(y)] = {"x_p": float(par["x_p"]), "x_e": float(par["x_e"]), "tau_days": int(par["tau"]),
                        "fit_days": int(par["n"])}
    OUT.write_text(yaml.safe_dump({
        "what": "SPAMS refitted on the design years leaving one year out (key = the year left out); X-107 amendment_3",
        "computed_on": str(date.today()), "settings": "insar_wetlands.spams.fit, freeze_v4_x106 mask and LAG",
        "frozen_reproduced": True, "design_years": list(F.DESIGN), "fits": fits}, sort_keys=False), encoding="utf-8")
    print(yaml.safe_dump(fits, sort_keys=False))


if __name__ == "__main__":
    main()
