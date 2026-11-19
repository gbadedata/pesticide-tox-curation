# Curation decisions

Each decision that changes or interprets a source value is listed here, is flagged on every
affected row through `issues`, and can be reversed by filtering on that code.

| # | Decision | Why | Flag on the row |
|---|---|---|---|
| 1 | Scope = substances in any OpenFoodTox dossier whose food domain is *pesticides*. Active substances (EFSA code `PPP`) and residues/metabolites (`PAR`) are kept and labelled `efsa_role`. | Mirrors how EFSA files the data; metabolites matter for residue definitions. | - |
| 2 | IUCLID `other:` + a free-text unit that is already a picklist unit is accepted as that unit. | 2,430 effect levels typed `mg/kg bw/day` as free text instead of choosing it. | `unit_picklist_bypassed` |
| 3 | `Âµ` is repaired to `µ`. | UTF-8 text decoded as Latin-1 upstream. | `unit_encoding_repaired` |
| 4 | ADI/AOEL recorded in `mg/kg bw` are read as `mg/kg bw/day`. | Both are daily values by definition; the `/day` was dropped. | `per_day_missing_from_unit` |
| 5 | LD50 recorded in `mg/kg bw/day` is read as `mg/kg bw`. | An LD50 is a single dose. | `per_day_on_single_dose` |
| 6 | LD50 recorded in `mg/kg` is read as `mg/kg bw`. | An LD50 is per kg body weight by definition. | `bw_assumed_for_ld50` |
| 7 | LD50 with a concentration unit (`mg/L`) is quarantined, not relabelled. | It is probably an LC50, but that is the data owner's call. | `ld50_with_concentration_unit` |
| 8 | `mg/kg`, `ppm` and per-person units on NOAELs and similar are not converted. | Converting needs food intake or body weight that the record does not give. | `unit_not_convertible` |
| 9 | A single value stored only in the upper-bound field (no qualifier) is used as the value. | IUCLID range fields misused for a point value. | `value_in_upper_field` |
| 10 | ARfD and AAOEL accept both `mg/kg bw` and `mg/kg bw/day`. | Both are single-day values; EFSA writes both forms. Canonical unit follows the record. | - |
| 11 | ChEBI roles are taken only from `has_role` (RO:0000087) asserted on the substance itself, not inherited from structural parents. | Inheritance over-generalises: glyphosate is_a *phosphonic acids*, which has the role *fungicide*. | - |
| 12 | ChEBI mapping uses the full InChIKey first, then a CAS number that maps to exactly one ChEBI entry. | A skeleton-only match can join a salt or isomer to the wrong entry. | `chebi_match` |
| 13 | AOP-Wiki matching falls back to the InChIKey skeleton (first block) only after full InChIKey and CAS fail. | Gains a few salts/esters; labelled so it can be excluded. | `match = inchikey_skeleton` |
| 14 | Common species names (rat, dog, ...) map to NCBI Taxonomy through a curated table whose targets must exist as NCBI scientific names. "monkey", "earthworm", "mammal" stay unmapped. | Regulatory records use common names; an unnamed species should not be guessed. | `taxon_match` |
| 15 | "Rounded" in the consistency check = implied value rounded to nearest at 1 or 2 significant figures, or up/down at 2. | Covers EFSA's rounding without excusing a twofold gap. | `status` |
| 16 | AOELs below NOAEL/UF with no oral absorption recorded are reported separately, not as unexplained. | EFSA corrects systemic AOELs for oral absorption, but the export rarely records the factor. | `status = absorption_not_recorded` |
