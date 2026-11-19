import pytest

from ptc import identifiers, units


@pytest.mark.parametrize("cas, ok", [
    ("1071-83-6", True),    # glyphosate
    ("2921-88-2", True),    # chlorpyrifos
    ("621-41-2", False),    # as recorded in OpenFoodTox for 3-acetamidophenol: bad check digit
    ("621-42-1", True),     # the valid CAS for 3-acetamidophenol
    ("not-a-cas", False),
    (None, False),
])
def test_cas_check_digit(cas, ok):
    assert identifiers.cas_is_valid(cas) is ok


def test_normalise_cas_strips_leading_zeros_and_text():
    assert identifiers.normalise_cas("0001071-83-6") == "1071-83-6"
    assert identifiers.normalise_cas("CAS 1071-83-6 (acid)") == "1071-83-6"
    assert identifiers.normalise_cas("n/a") is None


def test_other_picklist_with_standard_unit_is_flagged_and_resolved():
    unit, flags = units.clean_raw_unit("other:", "mg/kg bw/day")
    assert unit == "mg/kg bw/day"
    assert flags == ["unit_picklist_bypassed"]


def test_mojibake_is_repaired():
    unit, flags = units.clean_raw_unit("Âµg/kg bw/day", None)
    assert unit == "µg/kg bw/day"
    assert "unit_encoding_repaired" in flags


def test_weekly_dose_converted_to_daily():
    h = units.harmonise(7.0, "mg/kg bw/week", "TWI")
    assert h["value_canonical"] == pytest.approx(1.0)
    assert h["unit_canonical"] == "mg/kg bw/day"
    assert h["unit_ucum"] == "mg/kg/wk"


def test_microgram_converted():
    h = units.harmonise(30.0, "µg/kg bw/day", "ADI")
    assert h["value_canonical"] == pytest.approx(0.03)


def test_adi_without_per_day_is_harmonised_but_flagged():
    h = units.harmonise(0.1, "mg/kg bw", "ADI")
    assert h["unit_issue"] == "per_day_missing_from_unit"
    assert h["unit_canonical"] == "mg/kg bw/day"


def test_ld50_per_day_is_a_single_dose():
    h = units.harmonise(1700.0, "mg/kg bw/day", "LD50")
    assert h["unit_issue"] == "per_day_on_single_dose"
    assert h["unit_canonical"] == "mg/kg bw"
    assert h["value_canonical"] == 1700.0


def test_ld50_with_concentration_unit_is_not_converted():
    h = units.harmonise(2.7, "mg/L", "LD50")
    assert h["unit_issue"] == "ld50_with_concentration_unit"
    assert h["value_canonical"] is None


def test_dietary_units_are_never_converted():
    for u in ("mg/kg", "ppm", "mg/day"):
        h = units.harmonise(100.0, u, "NOAEL")
        assert h["value_canonical"] is None
        assert h["unit_issue"] == "unit_not_convertible"


def test_unknown_unit_is_not_guessed():
    assert units.harmonise(1.0, "mg/kg feed DM", "NOAEL")["unit_issue"] == "unit_unknown"
