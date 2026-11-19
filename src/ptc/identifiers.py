"""Chemical identifier checks and resolution against PubChem results.

Every substance ends up with exactly one ``id_status``:

resolved            one PubChem CID for the CAS number
resolved_parent     several CIDs, but they share one InChIKey first block (same skeleton,
                    e.g. stereo variants) - the lowest CID is kept and the rest recorded
ambiguous           several CIDs with different skeletons - not resolved, sent to review
not_found           PubChem has no record for the CAS number
invalid_cas         the CAS number fails its check digit
no_cas              the source gives no CAS number
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


def cas_is_valid(cas: str | None) -> bool:
    if not isinstance(cas, str) or not CAS_RE.match(cas):
        return False
    digits = cas.replace("-", "")
    body, check = digits[:-1], int(digits[-1])
    return sum(int(d) * w for w, d in enumerate(reversed(body), start=1)) % 10 == check


def normalise_cas(value) -> str | None:
    """Trim, drop leading zeros and reject anything that is not CAS-shaped."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip()
    if not s or s.lower() in ("nan", "none", "-", "n/a", "na"):
        return None
    m = re.search(r"\d{2,7}-\d{2}-\d", s)
    if not m:
        return None
    first, mid, last = m.group(0).split("-")
    return f"{int(first)}-{mid}-{last}"


def load_pubchem(raw: Path) -> tuple[dict[str, dict], pd.DataFrame, dict[int, list[str]]]:
    cas_rows = {}
    for line in (raw / "pubchem_cas_to_cid.jsonl").read_text().splitlines():
        r = json.loads(line)
        cas_rows[r["cas"]] = r
    props = pd.DataFrame([json.loads(line) for line in
                          (raw / "pubchem_properties.jsonl").read_text().splitlines()])
    chebi: dict[int, list[str]] = {}
    path = raw / "pubchem_cid_chebi.jsonl"
    if path.exists():
        for line in path.read_text().splitlines():
            r = json.loads(line)
            chebi[int(r["cid"])] = r["chebi"]
    return cas_rows, props, chebi


def resolve(cas: str | None, cas_rows: dict[str, dict], props: pd.DataFrame) -> dict:
    out = {"pubchem_cid": None, "inchikey": None, "id_status": None, "alt_cids": None}
    if cas is None:
        out["id_status"] = "no_cas"
        return out
    if not cas_is_valid(cas):
        out["id_status"] = "invalid_cas"
        return out
    row = cas_rows.get(cas)
    if row is None or not row["cids"]:
        out["id_status"] = "not_found"
        return out
    cids = sorted(set(int(c) for c in row["cids"]))
    keys = props.set_index("CID")["InChIKey"].reindex(cids)
    if len(cids) == 1:
        out.update(pubchem_cid=cids[0], inchikey=keys.iloc[0], id_status="resolved")
        return out
    skeletons = {k.split("-")[0] for k in keys.dropna()}
    if len(skeletons) == 1:
        out.update(pubchem_cid=cids[0], inchikey=keys.iloc[0], id_status="resolved_parent",
                   alt_cids=",".join(map(str, cids[1:])))
    else:
        out.update(id_status="ambiguous", alt_cids=",".join(map(str, cids)))
    return out
