"""Harmonise EFSA toxicological reference values (ADI, ARfD, AOEL, AAOEL, other).

Input: FLEX_SUM.ToxRefValues, one wide row per summary document with one block of
columns per descriptor. Output: one long row per (summary document, descriptor) with
a numeric value, UCUM unit, canonical value, uncertainty factor, the pivotal study
record it cites, the EFSA output it came from, and reason-coded issues.
"""
from __future__ import annotations

import pandas as pd

from . import units

P = "HumanHealthHazardCharacteristics."
BLOCKS = {
    "ADI": ("AcceptableDailyIntake", "Adi.lowerQualifier", "Adi.lowerValue",
            "Adi.upperQualifier", "Adi.upperValue", "Adi.Unit", None),
    "ARfD": ("AcuteReferenceDose", None, "Arfd.lowerValue",
             "Arfd.upperQualifier", "Arfd.upperValue", "Arfd.Unit", None),
    "AOEL": ("AcceptableOperatorExposureLevel", None, "Aoel.Value", None, None, "Aoel.Unit", None),
    "AAOEL": ("AcuteAcceptableOperatorExposureLevel", None, "Aaoel.Value", None, None,
              "Aaoel.Unit", None),
    "OTHER": ("OtherReferenceValues", "RefValue.lowerQualifier", "RefValue.lowerValue",
              "RefValue.upperQualifier", "RefValue.upperValue", "RefValue.Unit",
              "RefValue.Unit.Other"),
}
# Issues that make a row unusable as a numeric reference value.
HARD = {"value_missing", "value_not_numeric", "value_non_positive", "unit_missing",
        "unit_unknown", "unit_dimension_mismatch", "unit_other_without_text",
        "no_allocated_but_value", "ld50_with_concentration_unit"}


def _get(row, col):
    if col is None or col not in row.index:
        return None
    v = row[col]
    return None if v is None or (isinstance(v, float) and pd.isna(v)) else v


def _num(v):
    if v is None:
        return None, None
    try:
        return float(str(v).replace(",", ".").strip()), None
    except ValueError:
        return None, "value_not_numeric"


def harmonise(trv: pd.DataFrame, record_dossier: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    for _, r in trv.iterrows():
        for desc, (blk, lq, lv, uq, uv, unit_col, unit_other_col) in BLOCKS.items():
            b = P + blk + "."
            lower_raw = _get(r, b + lv)
            upper_raw = _get(r, b + uv) if uv else None
            no_alloc = str(_get(r, b + "NoAllocated")).lower() == "true"
            unit_raw = _get(r, b + unit_col)
            if lower_raw is None and upper_raw is None and not no_alloc:
                continue
            issues: list[str] = []
            descriptor = desc
            if desc == "OTHER":
                d = _get(r, b + "ReferenceValueDescriptor")
                descriptor = _get(r, b + "ReferenceValueDescriptor.Other") if d == "other:" else d
                descriptor = descriptor or "other"
            lower, e1 = _num(lower_raw)
            upper, e2 = _num(upper_raw)
            issues += [e for e in (e1, e2) if e]
            qualifier = _get(r, b + lq) if lq else None
            upper_q = _get(r, b + uq) if uq else None
            if lower is not None and upper is not None and lower != upper:
                value, value_field = lower, "range"
                issues.append("value_is_range")
            elif lower is not None:
                value, value_field = lower, "lower"
            elif upper is not None:
                value, value_field = upper, "upper"
                qualifier = upper_q or qualifier
                if not upper_q:
                    issues.append("value_in_upper_field")
            else:
                value, value_field = None, None
            if no_alloc and value is not None:
                issues.append("no_allocated_but_value")
            if value is None and not no_alloc and "value_not_numeric" not in issues:
                issues.append("value_missing")
            if value is not None and value <= 0:
                issues.append("value_non_positive")

            unit_clean, uflags = units.clean_raw_unit(unit_raw, _get(r, b + unit_other_col)
                                                      if unit_other_col else None)
            issues += uflags
            h = units.harmonise(value, unit_clean, descriptor) if not no_alloc else \
                {"unit_ucum": None, "dimension": None, "value_canonical": None,
                 "unit_canonical": None, "unit_issue": None}
            if h["unit_issue"] and not no_alloc:
                issues.append(h["unit_issue"])

            uf_raw = _get(r, b + "OverallUncertainty")
            uf, e3 = _num(uf_raw)
            if e3:
                issues.append("uf_not_numeric")
            if uf is None and value is not None and desc in ("ADI", "ARfD", "AOEL", "AAOEL"):
                issues.append("uf_missing")
            absorption, _ = _num(_get(r, b + "OralAbsorption"))
            rows.append({
                "record_uuid": r["Document UUID"],
                "substance_uuid": r["Parent UUID"],
                "descriptor": descriptor,
                "no_allocated": no_alloc,
                "qualifier": qualifier,
                "value": value,
                "value_upper": upper if value_field == "range" else None,
                "value_field": value_field,
                "unit_raw": unit_raw,
                "unit_clean": unit_clean,
                **h,
                "uncertainty_factor": uf,
                "oral_absorption_pct": absorption,
                "population": _get(r, b + "Population"),
                "critical_record_uuid": _get(r, b + "CriticalEndpoint"),
                "issues": ";".join(sorted(set(i for i in issues if i))),
            })
    out = pd.DataFrame(rows)
    out["quarantined"] = out["issues"].apply(
        lambda s: any(i in HARD for i in s.split(";")) if s else False)
    if record_dossier is not None:
        out = out.merge(record_dossier, on="record_uuid", how="left")
    return out
