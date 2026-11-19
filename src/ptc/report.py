"""Write reports/qa_report.md: aggregate curation results only (no EFSA data rows)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ISSUE_MEANING = {
    "unit_picklist_bypassed": "IUCLID 'other:' chosen, then a standard picklist unit typed as free text",
    "unit_free_text": "Unit given only as free text not in the picklist",
    "unit_only_in_free_text": "Unit field empty, unit found in the free-text field",
    "unit_encoding_repaired": "Unit text was UTF-8 read as Latin-1 ('Âµ' for 'µ'); repaired",
    "unit_other_without_text": "'other:' chosen with no free-text unit",
    "unit_missing": "No unit",
    "unit_unknown": "Unit string not in the curated unit table",
    "unit_not_convertible": "Unit needs extra data to convert (dietary ppm, mg/kg diet, per person)",
    "unit_dimension_mismatch": "Unit dimension impossible for the descriptor",
    "per_day_missing_from_unit": "ADI/AOEL given in mg/kg bw without '/day'; harmonised as per day",
    "per_day_on_single_dose": "Single-dose LD50 recorded per day; harmonised to mg/kg bw",
    "bw_assumed_for_ld50": "LD50 in 'mg/kg' read as mg/kg body weight",
    "ld50_with_concentration_unit": "LD50 with a concentration unit; probably an LC50",
    "value_in_upper_field": "Single value stored in the upper-bound field with no qualifier",
    "value_is_range": "Lower and upper values differ",
    "value_not_numeric": "Value is not a number",
    "value_non_positive": "Value is zero or negative",
    "uf_missing": "No overall uncertainty factor recorded",
    "uf_not_numeric": "Uncertainty factor is not a number",
    "no_allocated_but_value": "Marked 'not allocated' but a value is present",
    "species_too_vague": "Species recorded as 'mammal', 'other' or similar",
    "descriptor_missing": "Value with no descriptor (NOAEL, LD50, ...)",
}


def _issue_table(df: pd.DataFrame) -> str:
    s = df["issues"].fillna("").str.split(";").explode()
    s = _sorted_counts(s[s != ""])
    rows = ["| Issue | Rows | Meaning |", "|---|---:|---|"]
    rows += [f"| `{k}` | {v:,} | {ISSUE_MEANING.get(k, '')} |" for k, v in s.items()]
    return "\n".join(rows)


def _sorted_counts(series: pd.Series) -> pd.Series:
    """Counts by descending frequency, ties broken alphabetically (stable across pandas)."""
    vc = series.astype(str).value_counts()
    order = sorted(vc.index, key=lambda k: (-vc[k], k))
    return vc.reindex(order)


def _vc(series: pd.Series, header: tuple[str, str]) -> str:
    vc = _sorted_counts(series.fillna("(none)"))
    rows = [f"| {header[0]} | {header[1]} |", "|---|---:|"]
    rows += [f"| {k} | {v:,} |" for k, v in vc.items()]
    return "\n".join(rows)


def write(t: dict[str, pd.DataFrame], raw: Path, reports: Path) -> Path:
    subs, rv, el, cons, links = (t["substances"], t["reference_values"], t["effect_levels"],
                                 t["refvalue_consistency"], t["aop_links"])
    manifest = json.loads((raw / "manifest.json").read_text()) if (raw / "manifest.json").exists() else {"sources": {}}
    act = subs[subs["efsa_role"] == "active_substance"]
    num = rv[rv["descriptor"].isin(["ADI", "ARfD", "AOEL", "AAOEL"]) & rv["value_canonical"].notna()]
    ct = pd.crosstab(cons["descriptor"], cons["status"])
    order = [c for c in ["exact", "rounded", "absorption_not_recorded", "not_explained",
                         "not_checkable:no_usable_dose_in_cited_record",
                         "not_checkable:no_critical_record", "not_checkable:no_uncertainty_factor"]
             if c in ct.columns]
    ct = ct[order]
    ct_md = ["| Descriptor | " + " | ".join(order) + " | Total |",
             "|---|" + "---:|" * (len(order) + 1)]
    for d, row in ct.iterrows():
        ct_md.append(f"| {d} | " + " | ".join(f"{v:,}" for v in row) + f" | {row.sum():,} |")
    tot = ct.sum()
    ct_md.append("| **All** | " + " | ".join(f"**{v:,}**" for v in tot) + f" | **{tot.sum():,}** |")
    consistent = int(tot.get("exact", 0) + tot.get("rounded", 0))
    checkable = int(tot[[c for c in order if not c.startswith("not_checkable")]].sum())

    src_rows = ["| Source | File | Version / release | Licence | SHA-256 (first 12) |",
                "|---|---|---|---|---|"]
    for k, v in manifest["sources"].items():
        ver = v.get("doi") or v.get("release") or ""
        src_rows.append(f"| {k} | `{v['file']}` | {ver} | {v.get('licence', '')} | `{v['sha256'][:12]}` |")

    linked = links[links["aop_uuid"].notna()] if not links.empty else links
    endorsed = linked[linked["oecd_status"].fillna("").str.contains("Endorsed")] if not linked.empty else linked
    pubchem = ""
    if "cas_structure_check" in subs.columns:
        pubchem = ("\n### CAS number vs structure (PubChem)\n\nFor substances with a CAS number, "
                   "PubChem's structure for that CAS was compared with the curated InChIKey.\n\n"
                   + _vc(subs.loc[subs["cas"].notna(), "cas_structure_check"], ("Result", "Substances")))

    md = f"""# QA report: pesticide toxicology curation

Generated from the sources below. This report holds aggregate results only; the
curated tables stay local because OpenFoodTox is licensed CC BY-ND 4.0.

## Sources

{chr(10).join(src_rows)}

## Scope

Substances that appear in at least one OpenFoodTox dossier whose food domain is
*pesticides*: **{len(subs):,}** substances, of which **{len(act):,}** are active substances
(EFSA parameter code suffix `PPP`) and **{(subs['efsa_role'] == 'residue_or_metabolite').sum():,}**
are residues or metabolites (`PAR`).

## 1. Substance identity

{_vc(subs['structure_agreement'], ('InChI vs SMILES (RDKit InChIKeys)', 'Substances'))}

- CAS numbers: {subs['cas'].notna().sum():,} present, {int(subs['cas_valid'].eq(False).sum())} failing the check digit, {int(subs['cas_conflict'].sum())} conflicting between the two CAS fields.
- ChEBI: {subs['chebi_id'].notna().sum():,} substances mapped ({(subs['chebi_match'] == 'inchikey').sum():,} by InChIKey, {(subs['chebi_match'] == 'cas').sum():,} by CAS); {act['chebi_id'].notna().sum():,} of {len(act):,} active substances.
- Of the mapped active substances, {(act['chebi_id'].notna() & (act['chebi_classes'].fillna('') == '')).sum():,} have no agrochemical role asserted in ChEBI.
{pubchem}

### Active substances by ChEBI role

{_vc(act.loc[act['chebi_id'].notna(), 'chebi_classes'].replace('', '(no agrochemical role in ChEBI)'), ('ChEBI role(s)', 'Active substances'))}

## 2. Reference values (ADI, ARfD, AOEL, AAOEL and others)

{len(rv):,} reference-value rows; {len(num):,} numeric ADI/ARfD/AOEL/AAOEL values;
{int(rv['no_allocated'].sum())} marked "not allocated"; {int(rv['quarantined'].sum())} quarantined.

{_issue_table(rv)}

## 3. Study effect levels (NOAEL, LOAEL, BMDL, LD50, ...)

{len(el):,} effect-level rows from {el['record_uuid'].nunique():,} study records; {int(el['quarantined'].sum())} quarantined.

{_issue_table(el)}

Species mapping to NCBI Taxonomy:

{_vc(el['taxon_match'], ('Method', 'Rows'))}

## 4. Does each reference value follow from its cited study?

Each ADI/ARfD/AOEL/AAOEL cites a pivotal study record and an overall uncertainty factor
(UF). Expected value = point of departure in that record / UF (x oral absorption for
AOELs when recorded).

{chr(10).join(ct_md)}

**{consistent:,} of {checkable:,} checkable values ({consistent / checkable:.1%}) follow exactly or
after rounding.** Values not explained by the cited study, by pattern:

{_vc(cons.loc[cons['status'] == 'not_explained', 'pattern'], ('Pattern', 'Values'))}

`off_by_10^k` values differ from the implied value by an exact power of ten, which points
to a decimal slip or a UF recorded wrongly; `value_equals_point_of_departure` values equal
the study dose itself, as if no UF had been applied. Both are concrete queries for the
data owners.

## 5. Links to Adverse Outcome Pathways (AOP-Wiki)

{linked['substance_uuid'].nunique() if not linked.empty else 0} substances link to {linked['aop_uuid'].nunique() if not linked.empty else 0} AOPs
({endorsed['aop_uuid'].nunique() if not endorsed.empty else 0} OECD-endorsed), matched by InChIKey, then CAS, then InChIKey skeleton
(the last flagged as `inchikey_skeleton`, which can join a salt or stereoisomer to its parent).

{_vc(links['match'] if not links.empty else pd.Series(dtype=str), ('Match method', 'Links'))}

## Limitations

- Pesticide scope is defined by dossier food domain; a substance assessed only in other
  domains is out of scope even if it is used as a pesticide.
- ChEBI roles are used as asserted; missing roles are reported, not inferred.
- Effect levels cover human health records only; ecotoxicology records are not harmonised yet.
- Unit decisions that change a value's unit are flagged on every affected row and listed
  in `docs/decisions.md`.
"""
    path = reports / "qa_report.md"
    path.write_text(md)
    return path
