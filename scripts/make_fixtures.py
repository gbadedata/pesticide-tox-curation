#!/usr/bin/env python3
"""Cut the small, real test fixtures in tests/fixtures/ from the full raw sources.

Rows and elements are copied unmodified; only the selection is ours. Re-run after a
source update with:  python scripts/make_fixtures.py
"""
from __future__ import annotations

import copy
import gzip
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ptc import oft  # noqa: E402

RAW, FIX = ROOT / "data" / "raw", ROOT / "tests" / "fixtures"

# ToxRefValues summary documents chosen because each shows one behaviour under test.
TOXREF = {
    "46ac9b72-836e-468e-b1d5-9ad15fd9eb80": "ADI that follows exactly from its cited NOAEL / UF",
    "a76b5cf9-efc0-4200-a1ff-84670803c03e": "ADI that follows after rounding",
    "c8dafcc5-9d7c-4446-917d-7cd127e6a260": "AOEL below NOAEL/UF with no oral absorption recorded",
    "67ee9598-8e45-4cb5-845e-70d2bc6baa7b": "ADI off from NOAEL/UF by exactly a factor of ten",
    "153ab780-4f29-4561-aa0d-31bc450f02f0": "ADI stored in the upper-bound field",
    "8356198b-c7eb-4ddd-89ca-1a6570e893ed": "ADI in mg/kg bw without /day",
    "be26a34b-6958-4680-8f27-f2b1ce53b1ef": "unit with UTF-8/Latin-1 encoding damage",
}
CHEBI_COMPOUNDS = ["CHEBI:27744", "CHEBI:34631"]  # glyphosate, chlorpyrifos
CHEBI_EXTRA_ROLES = ["plant growth regulator", "herbicide safener", "agrochemical",
                     "molluscicide", "nematicide", "rodenticide"]
TAXIDS = {10116, 10090, 9615, 9986, 9606, 9031, 35525, 8022}
AOP_WIKI_ID = 3  # complex I inhibition -> parkinsonian motor deficits (rotenone)


def oft_slices() -> None:
    trv = oft.load("FLEX_SUM.ToxRefValues", RAW)
    t = trv[trv["Document UUID"].isin(TOXREF)]
    assert len(t) == len(TOXREF), "a chosen ToxRefValues document is missing"
    t.to_parquet(FIX / "oft_toxrefvalues.parquet", index=False)
    crit = set()
    for col in t.columns:
        if col.endswith(".CriticalEndpoint"):
            crit |= set(t[col].dropna())
    hh = oft.load("END_STUDY_REC.HumanHealth", RAW)
    hh[hh["Document UUID"].isin(crit)].to_parquet(FIX / "oft_humanhealth.parquet", index=False)


def chebi_slice() -> None:
    sys.path.insert(0, str(ROOT / "src"))
    from ptc.chebi import ChEBI
    c = ChEBI.from_obo(RAW / "chebi.obo.gz")
    keep: set[str] = set(CHEBI_COMPOUNDS)
    for cid in CHEBI_COMPOUNDS:
        keep.update(c.terms[cid].is_a)
    for t in list(keep):
        keep.update(c.terms[t].has_role)
    keep.update(filter(None, (c.role_id(n) for n in CHEBI_EXTRA_ROLES + ["pesticide"])))
    for t in list(keep):
        keep.update(c.ancestors(t))
    out, cur, take = [], [], False
    with gzip.open(RAW / "chebi.obo.gz", "rt", encoding="utf-8") as fh:
        header = True
        for line in fh:
            if line.startswith("["):
                if take:
                    out.extend(cur)
                header, cur, take = False, [line], False
                continue
            if header:
                out.append(line)
                continue
            cur.append(line)
            if line.startswith("id: "):
                take = line[4:].strip() in keep
        if take:
            out.extend(cur)
    (FIX / "chebi_slice.obo").write_text("".join(out), encoding="utf-8")


def taxonomy_slice() -> None:
    lines = []
    with gzip.open(RAW / "ncbi_taxonomy_names.dmp.gz", "rt", encoding="utf-8") as fh:
        for line in fh:
            if int(line.split("\t", 1)[0]) in TAXIDS:
                lines.append(line)
    (FIX / "ncbi_names_slice.dmp").write_text("".join(lines), encoding="utf-8")


def aop_slice() -> None:
    ns = "http://www.aopkb.org/aop-xml"
    ET.register_namespace("", ns)
    with gzip.open(next(RAW.glob("aop-wiki-xml*.gz"))) as fh:
        root = ET.parse(fh).getroot()

    def q(tag: str) -> str:
        return f"{{{ns}}}{tag}"

    vendor = root.find(q("vendor-specific"))
    aop_uuid = next(r.get("id") for r in vendor if r.tag == q("aop-reference")
                    and r.get("aop-wiki-id") == str(AOP_WIKI_ID))
    aop = next(e for e in root.findall(q("aop")) if e.get("id") == aop_uuid)
    stressors = {s.get("stressor-id") for s in aop.iter(q("aop-stressor"))}
    kes = {e.get("key-event-id") for e in aop.iter() if e.get("key-event-id")}
    stress_el = [s for s in root.findall(q("stressor")) if s.get("id") in stressors]
    chems = {c.get("chemical-id") for s in stress_el for c in s.iter(q("chemical-initiator"))}
    ke_el = [k for k in root.findall(q("key-event")) if k.get("id") in kes]
    bio = {v for k in ke_el for b in k.iter(q("biological-event")) for v in b.attrib.values()}
    keep_ids = {aop_uuid} | stressors | kes | chems | bio
    new = ET.Element(q("data"))
    for el in root:
        if el.tag == q("vendor-specific"):
            v = ET.SubElement(new, q("vendor-specific"))
            for r in el:
                if r.get("id") in keep_ids:
                    v.append(copy.deepcopy(r))
        elif el.get("id") in keep_ids:
            new.append(copy.deepcopy(el))
    ET.ElementTree(new).write(FIX / "aop_slice.xml", encoding="utf-8", xml_declaration=True)


if __name__ == "__main__":
    FIX.mkdir(parents=True, exist_ok=True)
    oft_slices()
    chebi_slice()
    taxonomy_slice()
    aop_slice()
    for p in sorted(FIX.iterdir()):
        print(f"{p.name}: {p.stat().st_size:,} bytes")
