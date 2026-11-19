# QA report: pesticide toxicology curation

Generated from the sources below. This report holds aggregate results only; the
curated tables stay local because OpenFoodTox is licensed CC BY-ND 4.0.

## Sources

| Source | File | Version / release | Licence | SHA-256 (first 12) |
|---|---|---|---|---|
| aopwiki | `aop-wiki-xml-2026-10-01.gz` | quarterly | see aopwiki.org terms | `317c591d2442` |
| chebi | `chebi.obo.gz` |  | CC BY 4.0 | `fdc1d2d5c4d3` |
| ncbi_taxonomy | `ncbi_taxonomy_names.dmp.gz` |  |  | `2d0ffdc3915f` |
| openfoodtox | `OFT3.0_export_repository.xlsx` | 10.5281/zenodo.19388272 | CC BY-ND 4.0 | `5181661e9216` |
| pubchem_cas_to_cid | `pubchem_cas_to_cid.jsonl` |  |  | `428d86347542` |
| pubchem_chebi | `pubchem_cid_chebi.jsonl` |  |  | `4345fdb601ef` |
| pubchem_properties | `pubchem_properties.jsonl` |  |  | `be4ebb7897f1` |

## Scope

Substances that appear in at least one OpenFoodTox dossier whose food domain is
*pesticides*: **1,513** substances, of which **410** are active substances
(EFSA parameter code suffix `PPP`) and **1,073**
are residues or metabolites (`PAR`).

## 1. Substance identity

| InChI vs SMILES (RDKit InChIKeys) | Substances |
|---|---:|
| full_match | 989 |
| no_structure | 376 |
| same_skeleton_different_stereo_or_charge | 94 |
| inchi_only | 25 |
| different_molecules | 13 |
| smiles_unparseable | 12 |
| smiles_only | 4 |

- CAS numbers: 738 present, 1 failing the check digit, 1 conflicting between the two CAS fields.
- ChEBI: 587 substances mapped (509 by InChIKey, 78 by CAS); 367 of 410 active substances.
- Of the mapped active substances, 112 have no agrochemical role asserted in ChEBI.

### CAS number vs structure (PubChem)

For substances with a CAS number, PubChem's structure for that CAS was compared with the curated InChIKey.

| Result | Substances |
|---|---:|
| agree | 613 |
| same_skeleton | 60 |
| not_found | 23 |
| no_structure_to_compare | 16 |
| same_parent_salt_or_mixture | 9 |
| metal_or_inorganic_representation | 8 |
| different_formula | 4 |
| same_formula_different_connectivity | 3 |
| ambiguous | 1 |
| invalid_cas | 1 |

### Active substances by ChEBI role

| ChEBI role(s) | Active substances |
|---|---:|
| (no agrochemical role in ChEBI) | 112 |
| herbicide | 80 |
| fungicide | 68 |
| insecticide | 26 |
| acaricide | 21 |
| acaricide;insecticide | 15 |
| plant growth regulator | 8 |
| herbicide;plant growth regulator | 6 |
| acaricide;insecticide;nematicide | 5 |
| agrochemical | 5 |
| nematicide | 5 |
| fungicide;nematicide | 3 |
| acaricide;insecticide;molluscicide | 2 |
| acaricide;nematicide | 2 |
| acaricide;fungicide | 1 |
| acaricide;insecticide;plant growth regulator | 1 |
| fungicide;herbicide | 1 |
| fungicide;herbicide;nematicide | 1 |
| fungicide;insecticide;nematicide | 1 |
| fungicide;insecticide;rodenticide | 1 |
| fungicide;plant growth regulator | 1 |
| insecticide;molluscicide | 1 |
| molluscicide | 1 |

## 2. Reference values (ADI, ARfD, AOEL, AAOEL and others)

2,613 reference-value rows; 2,197 numeric ADI/ARfD/AOEL/AAOEL values;
95 marked "not allocated"; 3 quarantined.

| Issue | Rows | Meaning |
|---|---:|---|
| `unit_not_convertible` | 186 | Unit needs extra data to convert (dietary ppm, mg/kg diet, per person) |
| `value_in_upper_field` | 159 | Single value stored in the upper-bound field with no qualifier |
| `uf_missing` | 104 | No overall uncertainty factor recorded |
| `unit_picklist_bypassed` | 47 | IUCLID 'other:' chosen, then a standard picklist unit typed as free text |
| `per_day_missing_from_unit` | 13 | ADI/AOEL given in mg/kg bw without '/day'; harmonised as per day |
| `uf_not_numeric` | 5 | Uncertainty factor is not a number |
| `unit_missing` | 2 | No unit |
| `unit_dimension_mismatch` | 1 | Unit dimension impossible for the descriptor |
| `unit_encoding_repaired` | 1 | Unit text was UTF-8 read as Latin-1 ('Âµ' for 'µ'); repaired |

## 3. Study effect levels (NOAEL, LOAEL, BMDL, LD50, ...)

5,248 effect-level rows from 3,010 study records; 13 quarantined.

| Issue | Rows | Meaning |
|---|---:|---|
| `unit_picklist_bypassed` | 2,430 | IUCLID 'other:' chosen, then a standard picklist unit typed as free text |
| `per_day_on_single_dose` | 462 | Single-dose LD50 recorded per day; harmonised to mg/kg bw |
| `species_too_vague` | 29 | Species recorded as 'mammal', 'other' or similar |
| `value_in_upper_field` | 25 | Single value stored in the upper-bound field with no qualifier |
| `bw_assumed_for_ld50` | 22 | LD50 in 'mg/kg' read as mg/kg body weight |
| `unit_not_convertible` | 18 | Unit needs extra data to convert (dietary ppm, mg/kg diet, per person) |
| `unit_missing` | 7 | No unit |
| `unit_other_without_text` | 7 | 'other:' chosen with no free-text unit |
| `ld50_with_concentration_unit` | 4 | LD50 with a concentration unit; probably an LC50 |
| `unit_dimension_mismatch` | 2 | Unit dimension impossible for the descriptor |
| `value_is_range` | 1 | Lower and upper values differ |

Species mapping to NCBI Taxonomy:

| Method | Rows |
|---|---:|
| curator_override | 5,209 |
| species_too_vague | 29 |
| missing | 7 |
| unique_other_name | 3 |

## 4. Does each reference value follow from its cited study?

Each ADI/ARfD/AOEL/AAOEL cites a pivotal study record and an overall uncertainty factor
(UF). Expected value = point of departure in that record / UF (x oral absorption for
AOELs when recorded).

| Descriptor | exact | rounded | absorption_not_recorded | not_explained | not_checkable:no_usable_dose_in_cited_record | not_checkable:no_critical_record | not_checkable:no_uncertainty_factor | Total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AAOEL | 81 | 18 | 20 | 3 | 29 | 2 | 3 | 156 |
| ADI | 473 | 173 | 0 | 102 | 18 | 52 | 14 | 832 |
| AOEL | 280 | 139 | 115 | 36 | 13 | 11 | 7 | 601 |
| ARfD | 407 | 58 | 0 | 46 | 55 | 33 | 9 | 608 |
| **All** | **1,241** | **388** | **135** | **187** | **115** | **98** | **33** | **2,197** |

**1,629 of 1,951 checkable values (83.5%) follow exactly or
after rounding.** Values not explained by the cited study, by pattern:

| Pattern | Values |
|---|---:|
| unexplained | 172 |
| off_by_10^-1 | 9 |
| off_by_10^1 | 5 |
| value_equals_point_of_departure | 1 |

`off_by_10^k` values differ from the implied value by an exact power of ten, which points
to a decimal slip or a UF recorded wrongly; `value_equals_point_of_departure` values equal
the study dose itself, as if no UF had been applied. Both are concrete queries for the
data owners.

## 5. Links to Adverse Outcome Pathways (AOP-Wiki)

43 substances link to 49 AOPs
(8 OECD-endorsed), matched by InChIKey, then CAS, then InChIKey skeleton
(the last flagged as `inchikey_skeleton`, which can join a salt or stereoisomer to its parent).

| Match method | Links |
|---|---:|
| inchikey | 71 |
| cas | 27 |
| inchikey_skeleton | 6 |

## Limitations

- Pesticide scope is defined by dossier food domain; a substance assessed only in other
  domains is out of scope even if it is used as a pesticide.
- ChEBI roles are used as asserted; missing roles are reported, not inferred.
- Effect levels cover human health records only; ecotoxicology records are not harmonised yet.
- Unit decisions that change a value's unit are flagged on every affected row and listed
  in `docs/decisions.md`.
