"""Unit harmonisation: raw OpenFoodTox unit strings -> UCUM code + canonical value.

Every raw unit seen in the pesticide subset is listed explicitly; an unknown unit is
never guessed and is reported as ``unit_unknown``.

Dimensions
----------
dose_rate         mass per kg body weight per day   canonical: mg/kg bw/day  (UCUM mg/kg/d)
dose              mass per kg body weight (single)  canonical: mg/kg bw      (UCUM mg/kg)
per_person_rate   mass per person per day            not converted (needs a body weight)
diet_or_dose      "mg/kg" with no "bw": dose or concentration in feed - not converted
diet_ppm          ppm in diet                        not converted (needs food intake)
water_conc        mass per litre                     canonical: mg/L
air_conc          mass per cubic metre               canonical: mg/m3
dimensionless
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Unit:
    ucum: str
    dimension: str
    canonical: str | None
    factor: float | None  # multiply raw value by this to get the canonical unit


UNITS: dict[str, Unit] = {
    "mg/kg bw/day": Unit("mg/kg/d", "dose_rate", "mg/kg bw/day", 1.0),
    "µg/kg bw/day": Unit("ug/kg/d", "dose_rate", "mg/kg bw/day", 1e-3),
    "ng/kg bw/day": Unit("ng/kg/d", "dose_rate", "mg/kg bw/day", 1e-6),
    "g/kg bw/day": Unit("g/kg/d", "dose_rate", "mg/kg bw/day", 1e3),
    "mg/kg bw/week": Unit("mg/kg/wk", "dose_rate", "mg/kg bw/day", 1 / 7),
    "µg/kg bw/week": Unit("ug/kg/wk", "dose_rate", "mg/kg bw/day", 1e-3 / 7),
    "mg/kg bw": Unit("mg/kg", "dose", "mg/kg bw", 1.0),
    "g/kg bw": Unit("g/kg", "dose", "mg/kg bw", 1e3),
    "mg/kg": Unit("mg/kg", "diet_or_dose", None, None),
    "g/kg": Unit("g/kg", "diet_or_dose", None, None),
    "ppm": Unit("[ppm]", "diet_ppm", None, None),
    "mg/day": Unit("mg/d", "per_person_rate", None, None),
    "mg/person/day": Unit("mg/d", "per_person_rate", None, None),
    "g/person/day": Unit("g/d", "per_person_rate", None, None),
    "mg/L": Unit("mg/L", "water_conc", "mg/L", 1.0),
    "µg/L": Unit("ug/L", "water_conc", "mg/L", 1e-3),
    "mg/m³": Unit("mg/m3", "air_conc", "mg/m3", 1.0),
    "dimensionless": Unit("1", "dimensionless", "1", 1.0),
}

# Text that was UTF-8 but decoded as Latin-1 somewhere upstream ("Âµ" for "µ").
MOJIBAKE = {"Âµ": "µ", "Ã‚Âµ": "µ", "mg/mÂ³": "mg/m³", "Â³": "³"}

# Which dimensions each reference descriptor may legitimately have.
EXPECTED = {
    "ADI": {"dose_rate"},
    "TDI": {"dose_rate"},
    "AOEL": {"dose_rate"},
    "ARfD": {"dose", "dose_rate"},   # a single-day value; EFSA writes both forms
    "AAOEL": {"dose", "dose_rate"},
    "NOAEL": {"dose_rate", "dose"},
    "LOAEL": {"dose_rate", "dose"},
    "NOEL": {"dose_rate", "dose"},
    "LOEL": {"dose_rate", "dose"},
    "BMDL10": {"dose_rate", "dose"},
    "BMDL05": {"dose_rate", "dose"},
    "BMD05": {"dose_rate", "dose"},
    "LD50": {"dose"},
}


def clean_raw_unit(unit, unit_other) -> tuple[str | None, list[str]]:
    """Resolve IUCLID's 'other:' picklist escape and repair encoding damage."""
    flags: list[str] = []
    u = None if unit is None or (isinstance(unit, float) and pd.isna(unit)) else str(unit).strip()
    o = None if unit_other is None or (isinstance(unit_other, float) and pd.isna(unit_other)) \
        else str(unit_other).strip()
    if u in (None, "") and o:
        u = o
        flags.append("unit_only_in_free_text")
    elif u == "other:":
        if o:
            u = o
            flags.append("unit_picklist_bypassed" if o in UNITS else "unit_free_text")
        else:
            return None, ["unit_other_without_text"]
    if u:
        for bad, good in MOJIBAKE.items():
            if bad in u:
                u = u.replace(bad, good)
                flags.append("unit_encoding_repaired")
    return (u or None), flags


def harmonise(value, unit: str | None, descriptor: str | None) -> dict:
    """Return UCUM code, dimension, canonical value/unit and any issue flags."""
    out = {"unit_ucum": None, "dimension": None, "value_canonical": None,
           "unit_canonical": None, "unit_issue": None}
    if unit is None:
        out["unit_issue"] = "unit_missing"
        return out
    spec = UNITS.get(unit)
    if spec is None:
        out["unit_issue"] = "unit_unknown"
        return out
    out.update(unit_ucum=spec.ucum, dimension=spec.dimension)
    factor, canonical = spec.factor, spec.canonical
    # Curation decisions below are flagged on the row and listed in docs/decisions.md.
    if descriptor == "LD50" and spec.dimension == "diet_or_dose":
        # an LD50 is per kg body weight by definition
        factor, canonical = (1.0 if spec.ucum == "mg/kg" else 1e3), "mg/kg bw"
        out["unit_issue"] = "bw_assumed_for_ld50"
    elif descriptor == "LD50" and spec.dimension == "dose_rate":
        # a single-dose LD50 recorded "per day": the number is a dose in mg/kg bw
        factor, canonical = spec.factor * (7 if spec.ucum.endswith("/wk") else 1), "mg/kg bw"
        out["unit_issue"] = "per_day_on_single_dose"
    elif descriptor == "LD50" and spec.dimension in ("water_conc", "air_conc"):
        out["unit_issue"] = "ld50_with_concentration_unit"  # probably an LC50
        return out
    elif spec.dimension in ("diet_or_dose", "diet_ppm", "per_person_rate"):
        # needs food intake or a body weight that the record does not give
        out["unit_issue"] = "unit_not_convertible"
        return out
    else:
        expected = EXPECTED.get(descriptor or "")
        if expected and spec.dimension not in expected:
            if descriptor in ("ADI", "TDI", "AOEL") and spec.dimension == "dose":
                out["unit_issue"] = "per_day_missing_from_unit"
                canonical = "mg/kg bw/day"
            else:
                out["unit_issue"] = "unit_dimension_mismatch"
                return out
    if factor is None or value is None or pd.isna(value):
        return out
    out["value_canonical"] = float(value) * factor
    out["unit_canonical"] = canonical
    return out
