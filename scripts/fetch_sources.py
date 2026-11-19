#!/usr/bin/env python3
"""Download every public source this project uses and record provenance.

Standard library only, so it runs on a fresh machine without installing anything:

    python3 scripts/fetch_sources.py --skip-pubchem          # files only (~5 min)
    python3 -m ptc.pipeline                                  # writes data/interim/pesticide_cas.txt
    python3 scripts/fetch_sources.py --pubchem-only --cas-list data/interim/pesticide_cas.txt
    python3 -m ptc.pipeline                                  # re-run with the PubChem cross-check

Re-running is safe: finished downloads are skipped (checksums are re-verified) and the
PubChem lookups resume from their cache.

Sources
-------
* EFSA OpenFoodTox 3.0 export (Zenodo, CC BY-ND 4.0) - not redistributed by this repo.
* AOP-Wiki XML export (latest quarterly release, falling back to the nightly file).
* ChEBI ontology (OBO) from EMBL-EBI.
* NCBI Taxonomy names (names.dmp from taxdump).
* PubChem PUG-REST: CAS number -> CID, then InChIKey / formula / SMILES / ChEBI IDs.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = RAW / "manifest.json"
UA = "pesticide-tox-curation/0.1 (github.com/gbadedata; research use)"

OFT_URL = ("https://zenodo.org/records/19388272/files/"
           "OFT3.0%20export%20repository.xlsx?download=1")
OFT_FILE = "OFT3.0_export_repository.xlsx"
OFT_MD5 = "445fc05a6a421634df822d14131a7d83"
OFT_DOI = "10.5281/zenodo.19388272"

AOP_DOWNLOADS_PAGE = "https://aopwiki.org/downloads"
AOP_NIGHTLY = "https://aopwiki.org/downloads/aop-wiki-xml.gz"

CHEBI_URLS = [
    "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz",
    "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo",
]
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CAS_RE = re.compile(r"\b(\d{2,7}-\d{2}-\d)\b")


# --------------------------------------------------------------------------- helpers
def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def http(url: str, data: bytes | None = None, timeout: int = 120, tries: int = 5) -> bytes:
    """GET (or POST when data is given) with retries on transient errors."""
    headers = {"User-Agent": UA}
    if data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    delay = 2.0
    for attempt in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code in (404, 400):
                raise
            if attempt == tries:
                raise
            log(f"  HTTP {e.code} on {url[:90]} - retry {attempt}/{tries - 1} in {delay:.0f}s")
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if attempt == tries:
                raise
            log(f"  {type(e).__name__} on {url[:90]} - retry {attempt}/{tries - 1} in {delay:.0f}s")
        time.sleep(delay)
        delay *= 2
    raise RuntimeError("unreachable")


def download(url: str, dest: Path) -> Path:
    tmp = dest.with_suffix(dest.suffix + ".part")
    headers = {"User-Agent": UA}
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as fh:
        total = int(resp.headers.get("Content-Length") or 0)
        done = 0
        last = time.time()
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            done += len(chunk)
            if time.time() - last > 10:
                pct = f" ({done / total:.0%})" if total else ""
                log(f"  {dest.name}: {done / 1e6:.1f} MB{pct}")
                last = time.time()
    tmp.replace(dest)
    return dest


def digests(path: Path) -> dict:
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
    return {"bytes": path.stat().st_size, "md5": md5.hexdigest(), "sha256": sha.hexdigest()}


def load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text())
    return {"sources": {}}


def save_manifest(m: dict) -> None:
    MANIFEST.write_text(json.dumps(m, indent=2, sort_keys=True))


def record(m: dict, key: str, path: Path, url: str, **extra) -> None:
    m["sources"][key] = {
        "file": path.name,
        "url": url,
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        **digests(path),
        **extra,
    }
    save_manifest(m)


# --------------------------------------------------------------------------- sources
def fetch_openfoodtox(m: dict) -> Path:
    dest = RAW / OFT_FILE
    if dest.exists() and digests(dest)["md5"] == OFT_MD5:
        log("OpenFoodTox already present, checksum OK")
    else:
        log("Downloading EFSA OpenFoodTox 3.0 export (22.6 MB)")
        download(OFT_URL, dest)
        got = digests(dest)["md5"]
        if got != OFT_MD5:
            sys.exit(f"OpenFoodTox MD5 mismatch: expected {OFT_MD5}, got {got}")
    record(m, "openfoodtox", dest, OFT_URL, doi=OFT_DOI, licence="CC BY-ND 4.0",
           expected_md5=OFT_MD5)
    return dest


def fetch_aopwiki(m: dict) -> Path:
    url, kind = AOP_NIGHTLY, "nightly"
    try:
        page = http(AOP_DOWNLOADS_PAGE).decode("utf-8", "replace")
        links = re.findall(r'href="([^"]*aop-wiki-xml-(\d{4}-\d{2}-\d{2})[^"]*\.gz)"', page)
        if links:
            href, _ = max(links, key=lambda t: t[1])
            url, kind = urllib.parse.urljoin(AOP_DOWNLOADS_PAGE, href), "quarterly"
    except Exception as e:  # noqa: BLE001 - fall back to nightly
        log(f"  could not read AOP-Wiki downloads page ({e}); using nightly export")
    dest = RAW / Path(urllib.parse.urlparse(url).path).name
    if not dest.exists():
        log(f"Downloading AOP-Wiki XML ({kind}): {url}")
        download(url, dest)
    else:
        log(f"AOP-Wiki XML already present: {dest.name}")
    record(m, "aopwiki", dest, url, release=kind, licence="see aopwiki.org terms")
    return dest


def fetch_chebi(m: dict) -> Path | None:
    for url in CHEBI_URLS:
        dest = RAW / Path(urllib.parse.urlparse(url).path).name
        if dest.exists():
            log(f"ChEBI already present: {dest.name}")
            record(m, "chebi", dest, url, licence="CC BY 4.0")
            return dest
        try:
            log(f"Downloading ChEBI ontology: {url}")
            download(url, dest)
            record(m, "chebi", dest, url, licence="CC BY 4.0")
            return dest
        except Exception as e:  # noqa: BLE001
            log(f"  failed ({e}); trying next ChEBI location")
    log("WARNING: ChEBI ontology not downloaded")
    return None


def fetch_taxonomy(m: dict) -> Path:
    dest = RAW / "ncbi_taxonomy_names.dmp.gz"
    if not dest.exists():
        log("Downloading NCBI Taxonomy taxdump (~65 MB) and keeping names.dmp only")
        tar_path = RAW / "taxdump.tar.gz"
        download(TAXDUMP_URL, tar_path)
        with tarfile.open(tar_path, "r:gz") as tar:
            member = tar.getmember("names.dmp")
            data = tar.extractfile(member).read()
        with gzip.open(dest, "wb") as fh:
            fh.write(data)
        tar_path.unlink()
    else:
        log("NCBI Taxonomy names already present")
    record(m, "ncbi_taxonomy", dest, TAXDUMP_URL, note="names.dmp extracted from taxdump.tar.gz")
    return dest


# --------------------------------------------------------------------------- PubChem
def cas_check_digit_ok(cas: str) -> bool:
    digits = cas.replace("-", "")
    body, check = digits[:-1], int(digits[-1])
    total = sum(int(d) * w for w, d in enumerate(reversed(body), start=1))
    return total % 10 == check


def cas_numbers_from_xlsx(path: Path) -> list[str]:
    """Every valid CAS number that appears as text anywhere in the workbook."""
    found: set[str] = set()
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            if name.endswith(".xml") and (name.startswith("xl/sharedStrings")
                                          or name.startswith("xl/worksheets/")):
                text = z.read(name).decode("utf-8", "replace")
                found.update(c for c in CAS_RE.findall(text) if cas_check_digit_ok(c))
    return sorted(found)


def pubchem(m: dict, oft: Path, cas_file: Path | None = None) -> None:
    if cas_file:
        cas_list = sorted({c.strip() for c in cas_file.read_text().split() if cas_check_digit_ok(c.strip())})
        log(f"PubChem: {len(cas_list):,} CAS numbers from {cas_file.name}")
    else:
        cas_list = cas_numbers_from_xlsx(oft)
        log(f"PubChem: {len(cas_list):,} valid CAS numbers found in OpenFoodTox")
    cache = RAW / "pubchem_cas_to_cid.jsonl"
    seen: dict[str, dict] = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            row = json.loads(line)
            seen[row["cas"]] = row
    todo = [c for c in cas_list if c not in seen]
    log(f"PubChem: {len(seen):,} cached, {len(todo):,} to look up (~{len(todo) / 4 / 60:.0f} min)")
    with open(cache, "a") as fh:
        for i, cas in enumerate(todo, 1):
            url = f"{PUG}/compound/name/{urllib.parse.quote(cas)}/cids/JSON"
            try:
                body = json.loads(http(url, timeout=60))
                row = {"cas": cas, "status": "found", "cids": body["IdentifierList"]["CID"]}
            except urllib.error.HTTPError as e:
                row = {"cas": cas, "status": "not_found" if e.code == 404 else f"http_{e.code}",
                       "cids": []}
            except Exception as e:  # noqa: BLE001
                row = {"cas": cas, "status": f"error:{type(e).__name__}", "cids": []}
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            seen[cas] = row
            if i % 200 == 0:
                log(f"  {i:,}/{len(todo):,} CAS looked up")
            time.sleep(0.25)  # stay under PubChem's 5 requests/second limit
    record(m, "pubchem_cas_to_cid", cache, f"{PUG}/compound/name/<CAS>/cids/JSON",
           n_cas=len(cas_list))

    cids = sorted({cid for r in seen.values() for cid in r["cids"]})
    log(f"PubChem: fetching properties for {len(cids):,} CIDs")
    props_path = RAW / "pubchem_properties.jsonl"
    wanted = ["InChIKey", "InChI", "MolecularFormula", "IUPACName", "SMILES", "ConnectivitySMILES"]
    fallback = ["InChIKey", "InChI", "MolecularFormula", "IUPACName", "IsomericSMILES",
                "CanonicalSMILES"]
    with open(props_path, "w") as fh:
        for start in range(0, len(cids), 300):
            chunk = cids[start:start + 300]
            data = urllib.parse.urlencode({"cid": ",".join(map(str, chunk))}).encode()
            try:
                body = json.loads(http(f"{PUG}/compound/cid/property/{','.join(wanted)}/JSON", data))
            except urllib.error.HTTPError:
                body = json.loads(http(f"{PUG}/compound/cid/property/{','.join(fallback)}/JSON",
                                       data))
            fh.writelines(json.dumps(p) + "\n" for p in body["PropertyTable"]["Properties"])
            time.sleep(0.3)
    record(m, "pubchem_properties", props_path, f"{PUG}/compound/cid/property/.../JSON")

    log("PubChem: extracting ChEBI IDs from synonyms")
    chebi_path = RAW / "pubchem_cid_chebi.jsonl"
    with open(chebi_path, "w") as fh:
        for start in range(0, len(cids), 200):
            chunk = cids[start:start + 200]
            data = urllib.parse.urlencode({"cid": ",".join(map(str, chunk))}).encode()
            body = json.loads(http(f"{PUG}/compound/cid/synonyms/JSON", data))
            for info in body["InformationList"]["Information"]:
                ids = sorted({s for s in info.get("Synonym", []) if re.fullmatch(r"CHEBI:\d+", s)})
                fh.write(json.dumps({"cid": info["CID"], "chebi": ids}) + "\n")
            time.sleep(0.3)
    record(m, "pubchem_chebi", chebi_path, f"{PUG}/compound/cid/synonyms/JSON",
           note="only CHEBI:<n> synonyms kept")


# --------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--skip-pubchem", action="store_true")
    ap.add_argument("--pubchem-only", action="store_true",
                    help="skip the file downloads and only run the PubChem step")
    ap.add_argument("--cas-list", type=Path,
                    help="look up only these CAS numbers (the pipeline writes "
                         "data/interim/pesticide_cas.txt)")
    args = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    m = load_manifest()
    if args.pubchem_only:
        oft = RAW / OFT_FILE
    else:
        oft = fetch_openfoodtox(m)
        fetch_aopwiki(m)
        fetch_chebi(m)
        fetch_taxonomy(m)
    if not args.skip_pubchem:
        pubchem(m, oft, args.cas_list)
    log(f"Done. Provenance written to {MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
