"""End-to-end build: raw sources -> curated tables, quarantine, data package, QA report.

    python -m ptc.pipeline            # writes data/curated/* and reports/qa_report.md
"""
from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from . import aopwiki, consistency, oft, refvalues, studies, substances
from .chebi import ChEBI
from .paths import CURATED, RAW, REPORTS
from .taxonomy import Taxonomy


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def record_provenance(raw: Path) -> pd.DataFrame:
    """Map each IUCLID document to the EFSA output(s) it was published in."""
    docs = oft.load("DOSSIER_DOCS", raw)
    dos = oft.dossiers(raw)
    m = docs.rename(columns={"DOSSIER UUID": "dossier_uuid", "DOCUMENT UUID": "record_uuid"})
    m = m[["record_uuid", "dossier_uuid"]].merge(dos, on="dossier_uuid", how="left")
    m["evaluation_date"] = pd.to_datetime(m["evaluation_date"], errors="coerce")
    m = m.sort_values("evaluation_date")
    g = m.groupby("record_uuid")
    return pd.DataFrame({
        "source_doi": g["doi"].agg(lambda s: ";".join(dict.fromkeys(s.dropna()))),
        "source_output_type": g["output_type"].agg(lambda s: ";".join(dict.fromkeys(s.dropna()))),
        "evaluation_date": g["evaluation_date"].max().dt.strftime("%Y-%m-%d"),
        "n_dossiers": g["dossier_uuid"].nunique(),
    }).reset_index()


def aop_links(subs: pd.DataFrame, aop: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Pesticide substance -> AOP-Wiki chemical -> stressor -> AOP -> MIE/AO key events."""
    chem = aop["chemicals"].copy()
    chem["key"] = chem["jchem_inchikey"].fillna(chem["indigo_inchikey"])
    links = []
    for _, s in subs.iterrows():
        k, cas = s.get("inchikey"), s.get("cas")
        hit, how = pd.DataFrame(), None
        if isinstance(k, str):
            hit, how = chem[chem["key"] == k], "inchikey"
        if hit.empty and isinstance(cas, str):
            hit, how = chem[chem["casrn"] == cas], "cas"
        if hit.empty and isinstance(k, str):
            hit, how = chem[chem["key"].str[:14] == k[:14]], "inchikey_skeleton"
        for _, c in hit.iterrows():
            links.append({"substance_uuid": s["substance_uuid"], "chemical_uuid": c["chemical_uuid"],
                          "aop_chemical_name": c["preferred_name"], "aop_casrn": c["casrn"],
                          "dsstox_id": c["dsstox_id"], "match": how})
    if not links:
        return pd.DataFrame()
    df = pd.DataFrame(links)
    df = (df.merge(aop["stressors"][["stressor_uuid", "stressor_name", "chemical_uuid"]],
                   on="chemical_uuid", how="left")
            .merge(aop["aop_stressor"][["aop_uuid", "stressor_uuid"]], on="stressor_uuid", how="left")
            .merge(aop["aops"], on="aop_uuid", how="left"))
    ke = aop["aop_ke"].merge(aop["key_events"], on="ke_uuid", how="left")
    for role in ("MIE", "AO"):
        sub = ke[ke["role"] == role].groupby("aop_uuid")["ke_title"].agg(
            lambda s: " | ".join(sorted(set(s.dropna()))))
        df[f"{role.lower()}_titles"] = df["aop_uuid"].map(sub)
    terms = ke.merge(aop["ke_terms"], on="ke_uuid", how="inner")
    terms = terms[terms["source"].isin(["GO", "UBERON", "CL", "HP", "MP", "CHEBI", "PR"])]
    tsum = terms.groupby("aop_uuid")["source_id"].agg(lambda s: ";".join(sorted(set(s))))
    df["ke_ontology_terms"] = df["aop_uuid"].map(tsum)
    return df


def run(raw: Path = RAW, curated: Path = CURATED, reports: Path = REPORTS) -> dict:
    curated.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    log("loading OpenFoodTox")
    pest = oft.pesticide_substance_uuids(raw)
    prov = record_provenance(raw)

    log("loading ChEBI and NCBI Taxonomy")
    chebi = ChEBI.from_obo(raw / "chebi.obo.gz")
    tax = Taxonomy.from_names_dmp(raw / "ncbi_taxonomy_names.dmp.gz")

    log("substance registry")
    subs_all = oft.substances(raw)
    subs = substances.build(subs_all[subs_all["substance_uuid"].isin(pest)], chebi, raw)

    cas_out = raw.parent / "interim" / "pesticide_cas.txt"
    cas_out.parent.mkdir(parents=True, exist_ok=True)
    cas_out.write_text("\n".join(sorted(subs["cas"].dropna().unique())) + "\n")

    log("reference values")
    trv = oft.load("FLEX_SUM.ToxRefValues", raw)
    rv = refvalues.harmonise(trv[trv["Parent UUID"].isin(pest)], prov)

    log("study effect levels")
    hh = oft.load("END_STUDY_REC.HumanHealth", raw)
    cited = set(rv["critical_record_uuid"].dropna())
    el = studies.harmonise(hh[hh["Parent UUID"].isin(pest) | hh["Document UUID"].isin(cited)], tax)
    el = el.merge(prov, on="record_uuid", how="left")

    log("consistency check")
    cons = consistency.check(rv, el)

    log("AOP-Wiki cross-map")
    aop = aopwiki.parse(next(raw.glob("aop-wiki-xml*.gz")))
    links = aop_links(subs, aop)

    quarantine = pd.concat([
        rv[rv["quarantined"]].assign(table="reference_values")[["table", "record_uuid", "substance_uuid", "descriptor", "issues"]],
        el[el["quarantined"]].assign(table="effect_levels")[["table", "record_uuid", "substance_uuid", "descriptor", "issues"]],
    ], ignore_index=True)

    tables = {"substances": subs, "reference_values": rv, "effect_levels": el,
              "refvalue_consistency": cons, "aop_links": links, "quarantine": quarantine}
    for name, df in tables.items():
        df.to_parquet(curated / f"{name}.parquet", index=False)
        df.to_csv(curated / f"{name}.csv", index=False)
    log(f"wrote {len(tables)} tables to {curated}")
    from . import package, report
    package.write(tables, raw, curated, curated.parents[1] / "schema" / "datapackage.json")
    report.write(tables, raw, reports)
    log(f"wrote {reports / 'qa_report.md'}")
    return tables


if __name__ == "__main__":
    run()
