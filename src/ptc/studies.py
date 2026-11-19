"""Harmonise human-health study effect levels (NOAEL, LOAEL, BMDL, LD50, ...).

IUCLID records hold effect levels in two structures that the export flattens:
``ResultsAndDiscussion.EffectLevels.Efflevel.*`` (repeated-dose style, NOAEL/LOAEL) and
``ResultsAndDiscussion.EffectLevels.*`` (acute style, LD50/LC50). Each export row is one
(study record, effect level) pair; both structures are read and the row is kept when
either holds a descriptor or a value.
"""
from __future__ import annotations

import re

import pandas as pd

from . import units
from .taxonomy import Taxonomy

E = "ResultsAndDiscussion.EffectLevels."
STRUCTS = {
    "efflevel": {"descriptor": E + "Efflevel.Endpoint", "descriptor_other": E + "Efflevel.Endpoint.Other",
                 "lq": E + "Efflevel.EffectLevel.lowerQualifier", "lv": E + "Efflevel.EffectLevel.lowerValue",
                 "uq": E + "Efflevel.EffectLevel.upperQualifier", "uv": E + "Efflevel.EffectLevel.upperValue",
                 "unit": E + "Efflevel.EffectLevel.Unit", "unit_other": E + "Efflevel.EffectLevel.Unit.Other",
                 "basis": E + "Efflevel.Basis", "sex": E + "Efflevel.Sex", "key": E + "Efflevel.KeyResult",
                 "uuid": E + "Efflevel.UUID"},
    "effectlevels": {"descriptor": E + "Endpoint", "descriptor_other": E + "Endpoint.Other",
                     "lq": E + "EffectLevel.lowerQualifier", "lv": E + "EffectLevel.lowerValue",
                     "uq": E + "EffectLevel.upperQualifier", "uv": E + "EffectLevel.upperValue",
                     "unit": E + "EffectLevel.Unit", "unit_other": E + "EffectLevel.Unit.Other",
                     "basis": E + "BasisForEffectLevel", "sex": E + "Sex", "key": None,
                     "uuid": E + "UUID"},
}
M = "MaterialsAndMethods."
HARD = {"value_not_numeric", "value_non_positive", "unit_missing", "unit_unknown",
        "unit_dimension_mismatch", "unit_other_without_text", "descriptor_missing",
        "ld50_with_concentration_unit"}
VAGUE_SPECIES = {"mammal", "other", "unspecified", "not specified", "rodents",
                 "mammalian cell line"}


def _g(r, col):
    if col is None or col not in r.index:
        return None
    v = r[col]
    return None if v is None or (isinstance(v, float) and pd.isna(v)) else str(v).strip()


def _num(v):
    if v is None:
        return None, None
    try:
        return float(v.replace(",", ".")), None
    except ValueError:
        return None, "value_not_numeric"


def species_text(r) -> str | None:
    s = _g(r, M + "TestAnimals.Species")
    if s in (None, "other:"):
        s = _g(r, M + "TestAnimals.Species.Other") or s
    return None if s in (None, "other:") else s


def harmonise(hh: pd.DataFrame, tax: Taxonomy | None) -> pd.DataFrame:
    rows = []
    for _, r in hh.iterrows():
        sp = species_text(r)
        taxid, tax_method = tax.resolve(sp) if tax else (None, "not_run")
        if sp and sp.lower() in VAGUE_SPECIES:
            taxid, tax_method = None, "species_too_vague"
        context = {
            "record_uuid": r["Document UUID"],
            "substance_uuid": r["Parent UUID"],
            "endpoint_category": _g(r, "AdministrativeData.Endpoint"),
            "study_type": _g(r, M + "StudyType"),
            "guideline": _g(r, M + "Guideline.Guideline.Other") or _g(r, M + "Guideline.Guideline"),
            "species_raw": sp,
            "ncbi_taxon": f"NCBITaxon:{taxid}" if taxid else None,
            "taxon_match": tax_method,
            "route": _g(r, M + "AdministrationExposure.RouteOfAdministration"),
            "duration": _g(r, M + "AdministrationExposure.DurationOfTreatmentExposure"),
        }
        for struct, c in STRUCTS.items():
            desc = _g(r, c["descriptor"])
            if desc in ("other:", "dose level:", "conc. level:"):
                desc = _g(r, c["descriptor_other"]) or desc.rstrip(":")
            lv_raw, uv_raw = _g(r, c["lv"]), _g(r, c["uv"])
            if desc is None and lv_raw is None and uv_raw is None:
                continue
            issues: list[str] = []
            lower, e1 = _num(lv_raw)
            upper, e2 = _num(uv_raw)
            issues += [e for e in (e1, e2) if e]
            if lower is not None:
                value, qualifier = lower, _g(r, c["lq"])
                if upper is not None and upper != lower:
                    issues.append("value_is_range")
            elif upper is not None:
                value, qualifier = upper, _g(r, c["uq"])
                if qualifier is None:
                    issues.append("value_in_upper_field")
            else:
                value, qualifier = None, None
                issues.append("value_missing")
            if desc is None:
                issues.append("descriptor_missing")
            if value is not None and value <= 0:
                issues.append("value_non_positive")
            unit_clean, uflags = units.clean_raw_unit(_g(r, c["unit"]), _g(r, c["unit_other"]))
            issues += uflags
            norm_desc = re.sub(r"\s+", "", desc).upper() if desc else None
            norm_desc = {"NOAEL": "NOAEL", "LOAEL": "LOAEL", "NOEL": "NOEL", "LOEL": "LOEL",
                         "BMDL10": "BMDL10", "BMDL05": "BMDL05", "BMD05": "BMD05",
                         "LD50": "LD50", "LC50": "LC50", "NOAEC": "NOAEC", "NOEC": "NOEC",
                         "LOAEC": "LOAEC"}.get(norm_desc or "", desc)
            h = units.harmonise(value, unit_clean, norm_desc) if value is not None else \
                {"unit_ucum": None, "dimension": None, "value_canonical": None,
                 "unit_canonical": None, "unit_issue": None}
            if h["unit_issue"]:
                issues.append(h["unit_issue"])
            if tax_method == "species_too_vague":
                issues.append("species_too_vague")
            rows.append({
                **context,
                "structure": struct,
                "effect_level_uuid": _g(r, c["uuid"]),
                "descriptor": norm_desc,
                "qualifier": qualifier,
                "value": value,
                "unit_raw": _g(r, c["unit"]),
                "unit_clean": unit_clean,
                **h,
                "basis": _g(r, c["basis"]),
                "sex": _g(r, c["sex"]),
                "key_result": (_g(r, c["key"]) or "").lower() == "true" if c["key"] else None,
                "issues": ";".join(sorted(set(issues))),
            })
    out = pd.DataFrame(rows)
    out["quarantined"] = out["issues"].apply(
        lambda s: any(i in HARD for i in s.split(";")) if s else False)
    # value_missing alone is common for records that only describe a study; keep them out
    # of the numeric table without calling them errors
    out["has_value"] = out["value"].notna()
    return out
