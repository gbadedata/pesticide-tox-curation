"""Substance registry: one row per pesticide-domain substance with curated identity.

Structure: InChIKeys are computed with RDKit from OpenFoodTox's InChI and SMILES; the
two must describe the same molecule (``structure_agreement``). The curated InChIKey
comes from the InChI when it parses, otherwise from the SMILES.

Identity cross-checks: CAS check digit, the two CAS fields in REF_SUB, ChEBI (by
InChIKey, then by CAS), and PubChem (CAS -> CID -> InChIKey) when its cache is present.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi as rd_inchi

from . import identifiers
from .chebi import ChEBI

RDLogger.DisableLog("rdApp.*")

EFSA_ROLE = {"PPP": "active_substance", "PAR": "residue_or_metabolite", "ADD": "additive",
             "VET": "veterinary", "CHE": "chemical"}


def _key_from_inchi(s) -> str | None:
    if not isinstance(s, str) or not s.startswith("InChI="):
        return None
    try:
        k = rd_inchi.InchiToInchiKey(s)
        return k or None
    except Exception:  # noqa: BLE001
        return None


def _key_from_smiles(s) -> str | None:
    if not isinstance(s, str) or not s.strip():
        return None
    m = Chem.MolFromSmiles(s.strip())
    if m is None:
        return None
    try:
        return Chem.MolToInchiKey(m) or None
    except Exception:  # noqa: BLE001
        return None


def _agreement(k_inchi, k_smiles, has_inchi, has_smiles) -> str:
    k_inchi = k_inchi if isinstance(k_inchi, str) else None
    k_smiles = k_smiles if isinstance(k_smiles, str) else None
    if has_inchi and k_inchi is None:
        return "inchi_unparseable"
    if has_smiles and k_smiles is None and not has_inchi:
        return "smiles_unparseable"
    if k_inchi and k_smiles:
        if k_inchi == k_smiles:
            return "full_match"
        if k_inchi[:14] == k_smiles[:14]:
            return "same_skeleton_different_stereo_or_charge"
        return "different_molecules"
    if k_inchi:
        return "inchi_only" if not has_smiles else "smiles_unparseable"
    if k_smiles:
        return "smiles_only"
    return "no_structure"


def build(subs: pd.DataFrame, chebi: ChEBI | None, raw: Path | None = None) -> pd.DataFrame:
    s = subs.copy()
    s["efsa_role"] = s["efsa_param_code"].str.extract(r"-([A-Z]+)$")[0].map(EFSA_ROLE)

    # CAS
    s["cas"] = s["cas_inventory"].map(identifiers.normalise_cas)
    s["cas_valid"] = s["cas"].map(lambda c: identifiers.cas_is_valid(c) if isinstance(c, str) else None)
    alt = s["cas_field"].map(identifiers.normalise_cas)
    s["cas_conflict"] = s["cas"].notna() & alt.notna() & (s["cas"] != alt)

    # Structure
    s["inchikey_from_inchi"] = s["inchi_oft"].map(_key_from_inchi)
    s["inchikey_from_smiles"] = s["smiles_oft"].map(_key_from_smiles)
    s["structure_agreement"] = [
        _agreement(a, b, isinstance(i, str), isinstance(m, str))
        for a, b, i, m in zip(s["inchikey_from_inchi"], s["inchikey_from_smiles"],
                              s["inchi_oft"], s["smiles_oft"], strict=True)]
    s["inchikey"] = s["inchikey_from_inchi"].fillna(s["inchikey_from_smiles"])

    # ChEBI
    if chebi is not None:
        def match(k, c):
            if isinstance(k, str) and k in chebi.by_inchikey:
                return chebi.by_inchikey[k][0], "inchikey"
            hits = chebi.by_cas.get(c, []) if isinstance(c, str) else []
            if len(hits) == 1:
                return hits[0], "cas"
            return None, ("cas_ambiguous" if hits else "no_match")
        pairs = [match(k, c) for k, c in zip(s["inchikey"], s["cas"], strict=True)]
        s["chebi_id"] = [p[0] for p in pairs]
        s["chebi_match"] = [p[1] for p in pairs]
        s["chebi_name"] = s["chebi_id"].map(lambda i: chebi.terms[i].name if isinstance(i, str) else None)
        s["chebi_classes"] = s["chebi_id"].map(
            lambda i: ";".join(chebi.pesticide_classes(i)) if isinstance(i, str) else None)

    # PubChem cross-check (optional)
    if raw is not None and (raw / "pubchem_cas_to_cid.jsonl").exists() \
            and (raw / "pubchem_properties.jsonl").exists():
        cas_rows, props, _ = identifiers.load_pubchem(raw)
        res = [identifiers.resolve(c, cas_rows, props) if isinstance(c, str) and c in cas_rows
               else {"pubchem_cid": None, "inchikey": None, "id_status": "not_looked_up",
                     "alt_cids": None}
               for c in s["cas"]]
        r = pd.DataFrame(res, index=s.index).add_prefix("pubchem_")
        s = pd.concat([s, r.rename(columns={"pubchem_pubchem_cid": "pubchem_cid"})], axis=1)

        def agree(row):
            a, b = row["inchikey"], row["pubchem_inchikey"]
            if not isinstance(b, str):
                return row["pubchem_id_status"]
            if not isinstance(a, str):
                return "no_structure_to_compare"
            if a == b:
                return "agree"
            return "same_skeleton" if a[:14] == b[:14] else "disagree"
        s["cas_structure_check"] = s.apply(agree, axis=1)
    return s
