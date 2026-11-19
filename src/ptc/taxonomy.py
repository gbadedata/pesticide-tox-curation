"""Map free-text test species to NCBI Taxonomy IDs using names.dmp.

Matching order:
1. curator overrides (common regulatory names -> a scientific name, documented in
   ``SPECIES_OVERRIDES``); the scientific name is then resolved through names.dmp so a
   wrong override fails loudly instead of producing a wrong ID;
2. exact (case-insensitive) scientific name;
3. a unique exact match on any NCBI name class (common name, synonym, ...).
Anything else is left unmapped with a reason, never guessed.
"""
from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path

# Regulatory toxicology uses common names for standard test species. Each entry is a
# curation decision; the right-hand side must exist as an NCBI scientific name.
SPECIES_OVERRIDES = {
    "rat": "Rattus norvegicus",
    "rats": "Rattus norvegicus",
    "mouse": "Mus musculus",
    "mice": "Mus musculus",
    "dog": "Canis lupus familiaris",
    "dogs": "Canis lupus familiaris",
    "rabbit": "Oryctolagus cuniculus",
    "rabbits": "Oryctolagus cuniculus",
    "hamster": "Mesocricetus auratus",
    "guinea pig": "Cavia porcellus",
    "monkey": None,  # genus/species not stated: leave unmapped on purpose
    "human": "Homo sapiens",
    "humans": "Homo sapiens",
    "chicken": "Gallus gallus",
    "hen": "Gallus gallus",
    "hens": "Gallus gallus",
    "pig": "Sus scrofa",
    "pigs": "Sus scrofa",
    "cattle": "Bos taurus",
    "cow": "Bos taurus",
    "goat": "Capra hircus",
    "sheep": "Ovis aries",
    "rainbow trout": "Oncorhynchus mykiss",
    "bluegill sunfish": "Lepomis macrochirus",
    "zebrafish": "Danio rerio",
    "honeybee": "Apis mellifera",
    "honey bee": "Apis mellifera",
    "bobwhite quail": "Colinus virginianus",
    "northern bobwhite": "Colinus virginianus",
    "mallard duck": "Anas platyrhynchos",
    "earthworm": None,  # several test species; leave unmapped unless named
}


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", str(s)).strip().lower()


class Taxonomy:
    def __init__(self, scientific: dict[str, int], any_name: dict[str, set[int]]):
        self.scientific = scientific
        self.any_name = any_name

    @classmethod
    def from_names_dmp(cls, path: Path) -> Taxonomy:
        opener = gzip.open if str(path).endswith(".gz") else open
        scientific: dict[str, int] = {}
        any_name: dict[str, set[int]] = defaultdict(set)
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 4:
                    continue
                taxid, name, _, klass = int(parts[0]), _norm(parts[1]), parts[2], parts[3]
                if klass == "scientific name":
                    scientific[name] = taxid
                any_name[name].add(taxid)
        return cls(scientific, dict(any_name))

    def resolve(self, text: str | None) -> tuple[int | None, str]:
        """Return (taxid, method) where method explains the match or why there is none."""
        if text is None or _norm(text) in ("", "nan", "none"):
            return None, "missing"
        key = _norm(text)
        if key in SPECIES_OVERRIDES:
            sci = SPECIES_OVERRIDES[key]
            if sci is None:
                return None, "unmapped_by_policy"
            taxid = self.scientific.get(_norm(sci))
            if taxid is None:
                raise ValueError(f"override {key!r} -> {sci!r} is not an NCBI scientific name")
            return taxid, "curator_override"
        if key in self.scientific:
            return self.scientific[key], "scientific_name"
        hits = self.any_name.get(key, set())
        if len(hits) == 1:
            return next(iter(hits)), "unique_other_name"
        if len(hits) > 1:
            return None, "ambiguous_name"
        return None, "no_match"
