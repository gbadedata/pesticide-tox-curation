# Test fixtures

Small slices of the real sources, cut by `scripts/make_fixtures.py`. Rows, OBO stanzas,
taxonomy lines and XML elements are copied unmodified; only the selection is ours.

| File | Source | Licence |
|---|---|---|
| `oft_toxrefvalues.parquet`, `oft_humanhealth.parquet` | EFSA OpenFoodTox 3.0 export, doi:10.5281/zenodo.19388272 (7 reference-value documents and the study records they cite, sheets converted to Parquet) | CC BY-ND 4.0, EFSA |
| `chebi_slice.obo` | ChEBI release 257 (glyphosate, chlorpyrifos, their parents and role lineages) | CC BY 4.0, EMBL-EBI |
| `ncbi_names_slice.dmp` | NCBI Taxonomy `names.dmp` (8 taxa) | Public domain, NCBI |
| `aop_slice.xml` | AOP-Wiki XML 2026-10-01 (AOP 3 and the stressors, chemicals and key events it references) | AOP-Wiki content licences as stated in each AOP |
