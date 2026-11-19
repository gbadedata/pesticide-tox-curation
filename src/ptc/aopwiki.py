"""Parse the AOP-Wiki XML export into tidy tables.

Tables produced
---------------
chemicals    one row per AOP-Wiki chemical (CAS, InChIKeys, DSSTox ID, preferred name)
stressors    stressor -> chemical links
aops         one row per AOP with its AOP-Wiki number, title and OECD status
aop_stressor stressor -> AOP links
aop_ke       AOP -> key event links, with the key event's role (MIE / KE / AO)
key_events   one row per key event with its level and any ontology terms
ke_terms     key event -> ontology term (GO, UBERON, CL, ...) from biological events
"""
from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _child(el: ET.Element, name: str) -> ET.Element | None:
    for c in el:
        if _local(c.tag) == name:
            return c
    return None


def _text(el: ET.Element, name: str) -> str | None:
    c = _child(el, name)
    if c is None or c.text is None:
        return None
    t = c.text.strip()
    return t or None


def _children(el: ET.Element, name: str) -> list[ET.Element]:
    c = _child(el, name)
    return [] if c is None else list(c)


def parse(path: Path) -> dict[str, pd.DataFrame]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        root = ET.parse(fh).getroot()

    chems, stressors, aops, aop_stress, aop_ke, kes, ke_terms = [], [], [], [], [], [], []
    bio: dict[str, dict] = {}  # biological object/process/action terms by uuid
    wiki_id: dict[tuple[str, str], int] = {}

    for el in root:
        tag = _local(el.tag)
        if tag in ("biological-object", "biological-process", "biological-action"):
            bio[el.get("id")] = {
                "kind": tag.replace("biological-", ""),
                "source": _text(el, "source"),
                "source_id": _text(el, "source-id"),
                "name": _text(el, "name"),
            }
        elif tag == "vendor-specific":
            for ref in el:
                rtag = _local(ref.tag).replace("-reference", "")
                if ref.get("aop-wiki-id"):
                    wiki_id[(rtag, ref.get("id"))] = int(ref.get("aop-wiki-id"))

    for el in root:
        tag = _local(el.tag)
        uid = el.get("id")
        if tag == "chemical":
            chems.append({
                "chemical_uuid": uid,
                "casrn": _text(el, "casrn"),
                "jchem_inchikey": _text(el, "jchem-inchi-key"),
                "indigo_inchikey": _text(el, "indigo-inchi-key"),
                "dsstox_id": _text(el, "dsstox-id"),
                "preferred_name": _text(el, "preferred-name"),
            })
        elif tag == "stressor":
            for ci in _children(el, "chemicals"):
                stressors.append({
                    "stressor_uuid": uid,
                    "stressor_name": _text(el, "name"),
                    "chemical_uuid": ci.get("chemical-id"),
                    "user_term": ci.get("user-term"),
                })
        elif tag == "key-event":
            kes.append({
                "ke_uuid": uid,
                "ke_wiki_id": wiki_id.get(("key-event", uid)),
                "ke_title": _text(el, "title"),
                "level": _text(el, "biological-organization-level"),
            })
            for be in _children(el, "biological-events"):
                for attr in ("object-id", "process-id", "action-id"):
                    term = bio.get(be.get(attr) or "")
                    if term and term["source_id"]:
                        ke_terms.append({"ke_uuid": uid, **term})
            for name in ("organ-term", "cell-term"):
                t = _child(el, name)
                if t is not None and _text(t, "source-id"):
                    ke_terms.append({"ke_uuid": uid, "kind": name.replace("-term", ""),
                                     "source": _text(t, "source"),
                                     "source_id": _text(t, "source-id"),
                                     "name": _text(t, "name")})
        elif tag == "aop":
            status = _child(el, "status")
            aops.append({
                "aop_uuid": uid,
                "aop_wiki_id": wiki_id.get(("aop", uid)),
                "aop_title": _text(el, "title"),
                "oecd_status": _text(status, "oecd-status") if status is not None else None,
                "wiki_status": _text(status, "wiki-status") if status is not None else None,
            })
            for s in _children(el, "aop-stressors"):
                aop_stress.append({"aop_uuid": uid, "stressor_uuid": s.get("stressor-id"),
                                   "evidence": _text(s, "evidence")})
            for c in el:
                ctag = _local(c.tag)
                if ctag == "molecular-initiating-event":
                    aop_ke.append({"aop_uuid": uid, "ke_uuid": c.get("key-event-id"), "role": "MIE"})
                elif ctag == "adverse-outcome":
                    aop_ke.append({"aop_uuid": uid, "ke_uuid": c.get("key-event-id"), "role": "AO"})
            for k in _children(el, "key-events"):
                aop_ke.append({"aop_uuid": uid, "ke_uuid": k.get("key-event-id") or k.get("id"),
                               "role": "KE"})

    aop_ke_df = pd.DataFrame(aop_ke, columns=["aop_uuid", "ke_uuid", "role"])
    # a key event listed as MIE/AO and again under key-events keeps the specific role
    order = {"MIE": 0, "AO": 1, "KE": 2}
    if not aop_ke_df.empty:
        aop_ke_df = (aop_ke_df.assign(_o=aop_ke_df["role"].map(order))
                     .sort_values("_o").drop_duplicates(["aop_uuid", "ke_uuid"])
                     .drop(columns="_o").reset_index(drop=True))
    return {
        "chemicals": pd.DataFrame(chems),
        "stressors": pd.DataFrame(stressors),
        "aops": pd.DataFrame(aops),
        "aop_stressor": pd.DataFrame(aop_stress),
        "aop_ke": aop_ke_df,
        "key_events": pd.DataFrame(kes),
        "ke_terms": pd.DataFrame(ke_terms).drop_duplicates() if ke_terms else pd.DataFrame(
            columns=["ke_uuid", "kind", "source", "source_id", "name"]),
    }
