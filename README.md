# pesticide-tox-curation

Curating EFSA's OpenFoodTox pesticide toxicology into FAIR, ontology-mapped tables, and
checking that each safety threshold actually follows from the study it cites.

All inputs are real public data: EFSA OpenFoodTox 3.0 (an IUCLID-structured export,
April 2026), AOP-Wiki (October 2026 release), ChEBI (release 257), NCBI Taxonomy and
PubChem. Nothing is synthetic, including the test fixtures.

## What it found

Scope: **1,513 substances** from **2,043 EFSA pesticide dossiers**: 410 active substances
and 1,073 residues or metabolites.

**Do the reference values follow from their cited studies?** Each ADI, ARfD, AOEL and
AAOEL in OpenFoodTox cites a pivotal study record and an uncertainty factor (UF), so the
value can be recomputed as *study dose / UF*.

| Result | Values |
|---|---:|
| Follow exactly or after rounding | **1,629 of 1,951 checkable (83.5%)** |
| AOELs explained by an oral-absorption correction that the record does not store | 135 |
| Not explained by the cited study | 187 |
| ...of which off by exactly 10x | 14 |
| ...of which equal to the study dose itself (no UF applied) | 1 |
| Not checkable (no cited study, no UF, or no usable dose in the cited study) | 246 |

The 15 patterned cases are specific queries for the data owner, each with its EFSA DOI
in the curated table. Examples: the ADI recorded for oxyfluorfen (0.003 mg/kg bw/day)
is one tenth of its cited NOAEL / UF (3 / 100), and the ARfD recorded for the metabolite
3,5,6-TCP equals its cited NOAEL (25 mg/kg bw). These have not yet been checked against
the published EFSA conclusions; they are leads, not verdicts.

**Chemical identity**

- 13 substances whose SMILES and InChI describe different molecules (InChIKeys computed
  with RDKit). For 3,5,6-TCP the InChI is trichloropyridinol (C5H2Cl3NO) but the SMILES
  is a dichloro, methyl-substituted molecule (C6H5Cl2NO). For a pyrazole metabolite
  named *...trifluoromethylsulfonyl...* the InChI matches the name, while the SMILES and
  the formula field both encode a sulfinyl group (one oxygen fewer).
- 94 more where the two agree on the skeleton but differ in stereochemistry or charge,
  and 12 SMILES that do not parse.
- One CAS number that fails its check digit (3-acetamidophenol, 621-41-2) and one
  substance with two different CAS numbers (triflumizole).
- 367 of 410 active substances map to ChEBI, but 112 of those have no agrochemical role
  asserted in ChEBI.
- PubChem cross-check (CAS number -> PubChem structure vs curated InChIKey): 613
  substances agree and 60 agree apart from stereochemistry. The 24 that don't are sorted
  into 9 salts or isomer mixtures, 8 inorganic or organometallic representations, 3 with
  the same formula but different connectivity (pyridalyl, spirodiclofen, tralkoxydim)
  and 4 with a different formula. For example, the CAS number recorded for the
  metabolite 2-methylpropyl methanesulfonate resolves in PubChem to a C11 terpene
  alcohol. The 23 CAS numbers PubChem doesn't know are mostly oils, extracts and
  mixtures with no single structure.

**Units**

- 2,430 effect levels where IUCLID's unit picklist was bypassed (`other:` chosen, then
  `mg/kg bw/day` typed by hand).
- 462 single-dose LD50s recorded "per day".
- 13 ADI/AOEL values missing "/day"; 159 reference values stored in the upper-bound
  field with no qualifier.
- One unit with encoding damage (`Âµg` for `µg`).

**Links to mechanisms**: 43 substances link to 49 Adverse Outcome Pathways in AOP-Wiki
(8 OECD-endorsed), carrying GO, UBERON, CL and other ontology terms from the AOPs' key
events.

The full breakdown is in [`reports/qa_report.md`](reports/qa_report.md).

## Pipeline

```
fetch_sources.py   OpenFoodTox 3.0 (Zenodo, MD5-checked) · AOP-Wiki XML · ChEBI OBO ·
                   NCBI Taxonomy names · PubChem (CAS -> CID -> InChIKey)
                   -> data/raw/ + manifest.json (URL, retrieval time, MD5, SHA-256, licence)

ptc.oft            IUCLID tables: DOSSIER -> SUBSTANCE -> REFERENCE_SUBSTANCE,
                   ToxRefValues and HumanHealth study records linked by document UUID
ptc.substances     CAS check digits · RDKit InChIKeys from InChI and SMILES ·
                   ChEBI by InChIKey/CAS · ChEBI agrochemical roles · PubChem cross-check
ptc.refvalues      ADI / ARfD / AOEL / AAOEL / other -> long table, UCUM units, reason codes
ptc.studies        NOAEL / LOAEL / BMDL / LD50 effect levels, species -> NCBI Taxonomy
ptc.consistency    reference value vs cited study dose / UF
ptc.aopwiki        substances -> AOP-Wiki chemicals -> stressors -> AOPs -> key events
ptc.package        Frictionless data package (schema/datapackage.json)
ptc.report         reports/qa_report.md (aggregate results only)
```

Every row keeps its IUCLID document UUID and the DOI of the EFSA output it came from.
Rows that cannot be used as numbers go to a quarantine table with reason codes, never
silently dropped.

## FAIR in practice

- **Findable**: every field is described in [`schema/datapackage.json`](schema/datapackage.json)
  (Frictionless Data Package), and every row carries persistent identifiers (IUCLID UUID,
  EFSA DOI, InChIKey, ChEBI, NCBI Taxonomy).
- **Accessible**: one command fetches every source from its public home, with checksums.
- **Interoperable**: units as UCUM codes, chemicals as InChIKey and ChEBI CURIEs, species
  as NCBITaxon CURIEs, mechanisms as AOP-Wiki IDs with GO and UBERON terms.
- **Reusable**: provenance and licence per source in `manifest.json`; every curation
  decision is listed in [`docs/decisions.md`](docs/decisions.md) and flagged on the rows
  it touches.

**Licensing.** OpenFoodTox is released under CC BY-ND 4.0, which does not allow sharing
modified versions. So this repository publishes the code, the schema and aggregate
results, and the curated tables are rebuilt locally from the original source. The test
fixtures are small, unmodified extracts with attribution
([`tests/fixtures/README.md`](tests/fixtures/README.md)).

## Run it

```bash
pip install -e ".[dev]"
python3 scripts/fetch_sources.py --skip-pubchem      # ~5 min, ~130 MB
python3 -m ptc.pipeline                              # ~1 min -> data/curated/, reports/
python3 scripts/fetch_sources.py --pubchem-only --cas-list data/interim/pesticide_cas.txt
python3 -m ptc.pipeline                              # adds the PubChem cross-check
```

Or `make fetch build pubchem build`. Outputs: `data/curated/*.parquet|csv` (substances,
reference_values, effect_levels, refvalue_consistency, aop_links, quarantine) and
`data/curated/datapackage.json`.

## Tests

`pytest` runs 29 tests on unmodified slices of the real sources: seven OpenFoodTox
reference-value documents chosen because each shows one behaviour (exact, rounded,
absorption, off-by-ten, value in the upper field, missing "/day", encoding damage),
the study records they cite, a ChEBI slice, an NCBI Taxonomy slice and the AOP-Wiki
entry for AOP 3. CI runs lint and tests on Python 3.10-3.12.

## Limitations

- Scope follows EFSA's dossier food domain; a substance assessed only in other domains
  is out of scope.
- Effect levels cover human-health study records. Ecotoxicology (aquatic and terrestrial
  records) is the next step.
- ChEBI roles are used as asserted; missing roles are reported, not inferred.
- The consistency check uses only the study EFSA cites. A value derived from another
  study, or from a UF the export records only in parts, shows as "not explained".
