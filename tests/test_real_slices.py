"""Tests on unmodified slices of the real sources (see tests/fixtures/README.md)."""
from pathlib import Path

import pandas as pd
import pytest

from ptc import aopwiki, consistency, refvalues, studies
from ptc.chebi import ChEBI
from ptc.taxonomy import SPECIES_OVERRIDES, Taxonomy

FIX = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="module")
def refvals():
    return refvalues.harmonise(pd.read_parquet(FIX / "oft_toxrefvalues.parquet"))


@pytest.fixture(scope="module")
def effect_levels():
    tax = Taxonomy.from_names_dmp(FIX / "ncbi_names_slice.dmp")
    return studies.harmonise(pd.read_parquet(FIX / "oft_humanhealth.parquet"), tax)


def _rv(refvals, record, descriptor="ADI"):
    r = refvals[(refvals["record_uuid"] == record) & (refvals["descriptor"] == descriptor)]
    assert len(r) == 1
    return r.iloc[0]


def test_adi_in_upper_field_is_read_and_flagged(refvals):
    r = _rv(refvals, "153ab780-4f29-4561-aa0d-31bc450f02f0")
    assert r["value_field"] == "upper"
    assert "value_in_upper_field" in r["issues"]
    assert r["value_canonical"] > 0


def test_adi_missing_per_day_is_flagged(refvals):
    r = _rv(refvals, "8356198b-c7eb-4ddd-89ca-1a6570e893ed")
    assert "per_day_missing_from_unit" in r["issues"]
    assert r["unit_canonical"] == "mg/kg bw/day"


def test_encoding_damage_repaired_on_real_row(refvals):
    rows = refvals[refvals["record_uuid"] == "be26a34b-6958-4680-8f27-f2b1ce53b1ef"]
    assert rows["issues"].str.contains("unit_encoding_repaired").any()
    assert not rows["unit_clean"].fillna("").str.contains("Â").any()


def test_consistency_statuses_on_real_records(refvals, effect_levels):
    c = consistency.check(refvals, effect_levels).set_index(["record_uuid", "descriptor"])
    assert c.loc[("46ac9b72-836e-468e-b1d5-9ad15fd9eb80", "ADI"), "status"] == "exact"
    assert c.loc[("a76b5cf9-efc0-4200-a1ff-84670803c03e", "ADI"), "status"] == "rounded"
    assert c.loc[("c8dafcc5-9d7c-4446-917d-7cd127e6a260", "AOEL"), "status"] == \
        "absorption_not_recorded"
    off = c.loc[("67ee9598-8e45-4cb5-845e-70d2bc6baa7b", "ADI")]
    assert off["status"] == "not_explained"
    assert off["pattern"] == "off_by_10^-1"


@pytest.mark.parametrize("x, sig, mode, expected", [
    (0.0145, 2, "half_up", 0.015),
    (0.016667, 2, "down", 0.016),
    (0.044228, 2, "up", 0.045),
    (0.125, 2, "half_up", 0.13),
])
def test_round_sig(x, sig, mode, expected):
    assert consistency.round_sig(x, sig, mode) == pytest.approx(expected)


def test_species_mapped_to_ncbi_taxonomy(effect_levels):
    rat = effect_levels[effect_levels["species_raw"].str.lower() == "rat"]
    assert not rat.empty
    assert set(rat["ncbi_taxon"]) == {"NCBITaxon:10116"}


def test_taxonomy_override_must_exist_in_ncbi():
    tax = Taxonomy.from_names_dmp(FIX / "ncbi_names_slice.dmp")
    assert tax.resolve("Rat") == (10116, "curator_override")
    assert tax.resolve("Daphnia magna") == (35525, "scientific_name")
    assert tax.resolve("monkey") == (None, "unmapped_by_policy")
    SPECIES_OVERRIDES["test-only"] = "Not a real species"
    try:
        with pytest.raises(ValueError):
            tax.resolve("test-only")
    finally:
        del SPECIES_OVERRIDES["test-only"]


def test_chebi_roles_are_direct_only():
    c = ChEBI.from_obo(FIX / "chebi_slice.obo")
    # a structural parent of glyphosate carries the fungicide role ...
    parents = c.terms["CHEBI:27744"].is_a
    assert any(c.pesticide_classes(p) == ["fungicide"] for p in parents)
    # ... but glyphosate itself is only asserted as a herbicide
    assert c.pesticide_classes("CHEBI:27744") == ["herbicide"]
    assert c.pesticide_classes("CHEBI:34631") == ["acaricide", "insecticide"]


def test_chebi_inchikey_index():
    c = ChEBI.from_obo(FIX / "chebi_slice.obo")
    assert c.by_inchikey["XDDAORKBJWWYJS-UHFFFAOYSA-N"] == ["CHEBI:27744"]


def test_aopwiki_slice_links_rotenone_to_aop3():
    t = aopwiki.parse(FIX / "aop_slice.xml")
    aop = t["aops"].iloc[0]
    assert aop["aop_wiki_id"] == 3
    assert aop["oecd_status"] == "WPHA/WNT Endorsed"
    assert "Rotenone" in set(t["chemicals"]["preferred_name"])
    roles = set(t["aop_ke"]["role"])
    assert {"MIE", "AO"} <= roles
    assert t["aop_ke"]["ke_uuid"].isin(t["key_events"]["ke_uuid"]).all()
