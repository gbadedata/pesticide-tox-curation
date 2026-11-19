"""Scientific consistency check: does each reference value follow from its cited study?

For ADI, ARfD, AOEL and AAOEL, EFSA records the pivotal study (CriticalEndpoint ->
END_STUDY_REC.HumanHealth) and the overall uncertainty factor (UF). The expected value is

    implied = point of departure (NOAEL / NOEL / BMDL / LOAEL in the cited record) / UF
              x oral absorption (AOEL/AAOEL, when recorded)

Each reference value gets exactly one status:

exact                      stated value within 1 % of implied
rounded                    stated value equals implied rounded to nearest at 1 or 2
                           significant figures, or up/down at 2 significant figures
absorption_not_recorded    AOEL/AAOEL below implied with no oral absorption recorded
                           (EFSA corrects systemic AOELs for oral absorption)
not_explained              none of the cited record's points of departure explain the value.
                           Causes include a wrong UF or decimal slip (see ``pattern``), a value
                           taken from a different study than the one cited, or a UF made of
                           parts the export does not record. It is a query, not a verdict.
not_checkable:<reason>     no critical record / no UF / cited record has no usable dose
"""
from __future__ import annotations

import math
from decimal import ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_UP, Decimal

import numpy as np
import pandas as pd

POD = ["NOAEL", "NOEL", "BMDL10", "BMDL05", "LOAEL", "NOAEC"]
CHECKED = ["ADI", "ARfD", "AOEL", "AAOEL"]


def round_sig(x: float, sig: int, mode: str = "half_up") -> float:
    """Round to ``sig`` significant figures: half_up (nearest), down or up."""
    if x == 0 or not math.isfinite(x):
        return x
    exp = math.floor(math.log10(abs(x))) - sig + 1
    q = Decimal(1).scaleb(exp)
    rounding = {"half_up": ROUND_HALF_UP, "down": ROUND_FLOOR, "up": ROUND_CEILING}[mode]
    return float(Decimal(repr(x)).quantize(q, rounding=rounding))


def _status(row) -> str:
    v, implied = row["value_canonical"], row["implied"]
    if abs(v / implied - 1) <= 0.01:
        return "exact"
    # to nearest at 1 or 2 significant figures; up or down only at 2 (a shift of at
    # most ~10 %), so a twofold gap is never excused as "rounding"
    for sig, mode in ((1, "half_up"), (2, "half_up"), (2, "down"), (2, "up")):
        r = round_sig(implied, sig, mode)
        if r and abs(v / r - 1) <= 0.01:
            return "rounded"
    if row["descriptor"] in ("AOEL", "AAOEL") and v < implied and pd.isna(row["oral_absorption_pct"]):
        return "absorption_not_recorded"
    return "not_explained"


def check(refvals: pd.DataFrame, effect_levels: pd.DataFrame) -> pd.DataFrame:
    rv = refvals[refvals["descriptor"].isin(CHECKED) & refvals["value_canonical"].notna()].copy()
    pod = effect_levels[effect_levels["descriptor"].isin(POD)
                        & effect_levels["value_canonical"].notna()
                        & ~effect_levels["quarantined"]]
    pod = pod[["record_uuid", "descriptor", "value_canonical", "species_raw"]].rename(columns={
        "record_uuid": "critical_record_uuid", "descriptor": "pod_descriptor",
        "value_canonical": "pod_value", "species_raw": "pod_species"})

    base = rv[["record_uuid", "substance_uuid", "descriptor", "value_canonical",
               "uncertainty_factor", "oral_absorption_pct", "critical_record_uuid"]]
    m = base.merge(pod, on="critical_record_uuid", how="left")
    for col in ("value_canonical", "uncertainty_factor", "oral_absorption_pct", "pod_value"):
        m[col] = pd.to_numeric(m[col], errors="coerce").astype(float)
    ok = m["pod_value"].notna() & m["uncertainty_factor"].notna() & (m["uncertainty_factor"] > 0)
    absorb = np.where(m["descriptor"].isin(["AOEL", "AAOEL"]) & m["oral_absorption_pct"].notna(),
                      m["oral_absorption_pct"] / 100, 1.0)
    m["implied"] = np.where(ok, m["pod_value"] / m["uncertainty_factor"] * absorb, np.nan)
    m["ratio"] = m["value_canonical"] / m["implied"]
    m["status"] = None
    m.loc[ok, "status"] = m[ok].apply(_status, axis=1)
    rank = {"exact": 0, "rounded": 1, "absorption_not_recorded": 2, "not_explained": 3}
    m["_rank"] = m["status"].map(rank).fillna(9)
    m["_dist"] = np.abs(np.log(m["ratio"].where(m["ratio"] > 0)))
    best = (m.sort_values(["_rank", "_dist"])
            .drop_duplicates(["record_uuid", "descriptor"])
            .drop(columns=["_rank", "_dist"]))

    def reason(r):
        if pd.notna(r["status"]):
            return r["status"]
        if pd.isna(r["critical_record_uuid"]):
            return "not_checkable:no_critical_record"
        if pd.isna(r["uncertainty_factor"]):
            return "not_checkable:no_uncertainty_factor"
        return "not_checkable:no_usable_dose_in_cited_record"

    best["status"] = best.apply(reason, axis=1)

    def pattern(r):
        if r["status"] != "not_explained":
            return None
        ratio = r["ratio"]
        if abs(r["value_canonical"] / r["pod_value"] - 1) <= 0.01:
            return "value_equals_point_of_departure"
        k = round(math.log10(ratio))
        if k != 0 and abs(ratio / 10 ** k - 1) <= 0.02:
            return f"off_by_10^{k}"
        return "unexplained"
    best["pattern"] = best.apply(pattern, axis=1)
    return best.reset_index(drop=True)
