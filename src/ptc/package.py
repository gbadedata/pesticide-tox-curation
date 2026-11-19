"""Frictionless Data Package descriptor for the curated tables (the FAIR metadata layer).

The descriptor is written next to the local data (data/curated/datapackage.json) and a
copy without row counts is kept in the repository (schema/datapackage.json), so the
structure, meaning, standards and licences of every field are findable even though the
OpenFoodTox-derived rows themselves cannot be redistributed (CC BY-ND 4.0).
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

TABLE_DOC = {
    "substances": "One row per pesticide-domain substance: identifiers, structure checks, ChEBI mapping, PubChem cross-check.",
    "reference_values": "Health-based reference values (ADI, ARfD, AOEL, AAOEL, other) in long format with harmonised units.",
    "effect_levels": "Human-health study effect levels (NOAEL, LOAEL, BMDL, LD50, ...) with study context and harmonised units.",
    "refvalue_consistency": "Whether each ADI/ARfD/AOEL/AAOEL follows from the point of departure and uncertainty factor of the study it cites.",
    "aop_links": "Links from substances to AOP-Wiki chemicals, stressors, AOPs and key-event ontology terms.",
    "quarantine": "Rows excluded from numeric use, with reason codes.",
}

FIELD_DOC = {
    "substance_uuid": ("IUCLID SUBSTANCE document UUID in OpenFoodTox 3.0", None),
    "record_uuid": ("IUCLID document UUID of the source record", None),
    "inchikey": ("Curated InChIKey, computed with RDKit from the InChI (else SMILES)", "https://www.inchi-trust.org/"),
    "chebi_id": ("ChEBI identifier (CURIE)", "http://purl.obolibrary.org/obo/chebi.owl"),
    "chebi_classes": ("Agrochemical classes asserted in ChEBI via has_role (RO:0000087)", "http://purl.obolibrary.org/obo/RO_0000087"),
    "cas": ("CAS Registry Number, normalised", None),
    "cas_valid": ("CAS check digit passes", None),
    "structure_agreement": ("Whether OpenFoodTox InChI and SMILES describe the same molecule", None),
    "pubchem_cid": ("PubChem compound ID found for the CAS number", "https://pubchem.ncbi.nlm.nih.gov/"),
    "cas_structure_check": ("PubChem structure for the CAS vs curated InChIKey", None),
    "efsa_role": ("active_substance (EFSA code PPP) or residue_or_metabolite (PAR), etc.", None),
    "descriptor": ("Reference value or dose descriptor (ADI, ARfD, AOEL, NOAEL, LD50, ...)", None),
    "value": ("Numeric value as recorded", None),
    "unit_raw": ("Unit as recorded in IUCLID", None),
    "unit_ucum": ("Unit as a UCUM code", "https://ucum.org/"),
    "value_canonical": ("Value converted to the canonical unit", None),
    "unit_canonical": ("Canonical unit: mg/kg bw/day, mg/kg bw, mg/L or mg/m3", None),
    "uncertainty_factor": ("Overall uncertainty factor applied by EFSA", None),
    "critical_record_uuid": ("Study record EFSA cites as the critical endpoint", None),
    "ncbi_taxon": ("Test species as an NCBI Taxonomy CURIE", "http://purl.obolibrary.org/obo/ncbitaxon.owl"),
    "taxon_match": ("How the species was mapped (curator_override, scientific_name, ...)", None),
    "issues": ("Semicolon-separated reason codes (see reports/qa_report.md)", None),
    "quarantined": ("True when an issue makes the row unusable as a number", None),
    "status": ("Consistency status: exact, rounded, absorption_not_recorded, not_explained, not_checkable:*", None),
    "pattern": ("For not_explained values: off_by_10^k, value_equals_point_of_departure or unexplained", None),
    "source_doi": ("DOI of the EFSA output(s) the record was published in", "https://doi.org"),
    "aop_wiki_id": ("AOP number in AOP-Wiki", "https://aopwiki.org"),
    "ke_ontology_terms": ("GO, UBERON, CL, HP, MP, CHEBI, PR terms on the AOP's key events", None),
}

LICENCES = [
    {"name": "CC-BY-ND-4.0", "title": "EFSA OpenFoodTox (derived rows are not redistributed)",
     "path": "https://creativecommons.org/licenses/by-nd/4.0/"},
    {"name": "CC-BY-4.0", "title": "ChEBI", "path": "https://creativecommons.org/licenses/by/4.0/"},
]


def _type(dtype) -> str:
    if pd.api.types.is_bool_dtype(dtype):
        return "boolean"
    if pd.api.types.is_integer_dtype(dtype):
        return "integer"
    if pd.api.types.is_float_dtype(dtype):
        return "number"
    return "string"


def descriptor(tables: dict[str, pd.DataFrame], manifest: dict, with_counts: bool) -> dict:
    resources = []
    for name, df in tables.items():
        fields = []
        for col in df.columns:
            f = {"name": col, "type": _type(df[col].dtype)}
            if col in FIELD_DOC:
                f["description"], ref = FIELD_DOC[col]
                if ref:
                    f["rdfType"] = ref
            fields.append(f)
        res = {"name": name, "path": f"{name}.csv", "format": "csv",
               "description": TABLE_DOC.get(name, ""), "schema": {"fields": fields}}
        if with_counts:
            res["rows"] = len(df)
        resources.append(res)
    sources = [{"title": k, "path": v["url"], "sha256": v["sha256"],
                **({"doi": v["doi"]} if "doi" in v else {})}
               for k, v in manifest.get("sources", {}).items()]
    return {"name": "pesticide-tox-curation", "title": "Curated EFSA pesticide toxicology",
            "profile": "data-package", "licenses": LICENCES, "sources": sources,
            "resources": resources}


def write(tables: dict[str, pd.DataFrame], raw: Path, curated: Path, repo_schema: Path) -> None:
    mpath = raw / "manifest.json"
    manifest = json.loads(mpath.read_text()) if mpath.exists() else {}
    (curated / "datapackage.json").write_text(json.dumps(descriptor(tables, manifest, True), indent=2))
    repo_schema.parent.mkdir(parents=True, exist_ok=True)
    repo_schema.write_text(json.dumps(descriptor(tables, manifest, False), indent=2))
