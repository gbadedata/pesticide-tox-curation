"""Load the EFSA OpenFoodTox 3.0 export (an IUCLID-structured workbook).

How the IUCLID tables link (verified on the 2026 export):

    DOSSIER --(DOSSIER_DOCS, type SUBSTANCE)--> SUB --(ReferenceSubstance)--> REF_SUB
    SUB <--(Parent UUID)-- FLEX_SUM.ToxRefValues, END_STUDY_REC.* , END_SUM
    ToxRefValues.*.CriticalEndpoint --> END_STUDY_REC.HumanHealth.Document UUID

Sheets are cached as Parquet on first load (the workbook takes ~30 s to parse).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .paths import INTERIM, RAW

OFT_FILE = "OFT3.0_export_repository.xlsx"
SHEETS_USED = [
    "DOSSIER", "DOSSIER_DOCS", "SUB", "REF_SUB", "FLEX_SUM.ToxRefValues",
    "END_STUDY_REC.HumanHealth",
]


def load(sheet: str, raw: Path = RAW, interim: Path = INTERIM) -> pd.DataFrame:
    cache = interim / f"{sheet}.parquet"
    if not cache.exists():
        interim.mkdir(parents=True, exist_ok=True)
        book = pd.read_excel(raw / OFT_FILE, sheet_name=None, dtype=str, engine="openpyxl")
        for name, df in book.items():
            df.drop(columns=[df.columns[0]]).to_parquet(interim / f"{name}.parquet", index=False)
    return pd.read_parquet(cache)


def dossiers(raw: Path = RAW, interim: Path = INTERIM) -> pd.DataFrame:
    d = load("DOSSIER", raw, interim)
    return pd.DataFrame({
        "dossier_uuid": d["Document UUID"],
        "food_domain": d["Domain.FoodDomain"],
        "regulation": d["Domain.Regulation"],
        "efsa_question": d["DataSource.EFSAQuestionNumber"],
        "output_title": d["LiteratureReference.EFSAOutputTitle"],
        "output_type": d["LiteratureReference.Type"],
        "evaluation_date": d["LiteratureReference.DateOfEvaluation"],
        "doi": d["LiteratureReference.LinkToPersistentIdentifier"].str.replace(
            r"^doi:\s*", "", regex=True),
    })


def substance_dossier_links(raw: Path = RAW, interim: Path = INTERIM) -> pd.DataFrame:
    docs = load("DOSSIER_DOCS", raw, interim)
    s = docs[docs["DOCUMENT TYPE"] == "SUBSTANCE"]
    return pd.DataFrame({"dossier_uuid": s["DOSSIER UUID"], "substance_uuid": s["DOCUMENT UUID"]})


def pesticide_substance_uuids(raw: Path = RAW, interim: Path = INTERIM) -> set[str]:
    """Substances that appear in at least one dossier whose food domain is pesticides."""
    d = dossiers(raw, interim)
    pest = set(d.loc[d["food_domain"] == "pesticides", "dossier_uuid"])
    links = substance_dossier_links(raw, interim)
    return set(links.loc[links["dossier_uuid"].isin(pest), "substance_uuid"])


def substances(raw: Path = RAW, interim: Path = INTERIM) -> pd.DataFrame:
    """One row per IUCLID SUBSTANCE joined to its REFERENCE_SUBSTANCE identifiers."""
    sub = load("SUB", raw, interim)
    ref = load("REF_SUB", raw, interim)
    s = pd.DataFrame({
        "substance_uuid": sub["Document UUID"],
        "chemical_name": sub["ChemicalName"],
        "composition": sub["TypeOfSubstance.Composition"],
        "ref_uuid": sub["ReferenceSubstance.ReferenceSubstance"],
    })
    r = pd.DataFrame({
        "ref_uuid": ref["Document UUID"],
        "reference_name": ref["ReferenceSubstanceName"],
        "cas_inventory": ref["Inventory.CASNumber"],
        "cas_field": ref["CAS number"],
        "ec_number": ref["EC number"],
        "iupac_name": ref["IupacName"],
        "inchikey_oft": ref["MolecularStructuralInfo.InChIKey"],
        "inchi_oft": ref["MolecularStructuralInfo.InChl"],
        "formula_oft": ref["MolecularStructuralInfo.MolecularFormula"],
        "smiles_oft": ref["MolecularStructuralInfo.SmilesNotation"],
        "pubchem_cid_oft": ref["PUBCHEM CID"],
        "dtxsid": ref["US EPA DTXSID"],
        "efsa_param_code": ref["EFSA PARAM CODE"],
    })
    return s.merge(r, on="ref_uuid", how="left", validate="many_to_one")
