"""ChEBI ontology: parse the OBO file and classify substances by pesticide role.

Role terms are looked up by name at run time (not hard-coded IDs). Only a substance's
own ``has_role`` (RO:0000087) assertions are used: inheriting roles from structural
parents over-generalises (glyphosate is_a phosphonic acid, and phosphonic acid has the
role fungicide). A role counts as agrochemical when ``pesticide``, ``agrochemical`` or
``plant growth regulator`` is in the role's own is_a lineage.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

PESTICIDE_CLASSES = [
    "herbicide", "fungicide", "insecticide", "acaricide", "nematicide",
    "molluscicide", "rodenticide", "plant growth regulator", "herbicide safener",
]


@dataclass
class Term:
    id: str
    name: str = ""
    is_a: list[str] = field(default_factory=list)
    has_role: list[str] = field(default_factory=list)
    obsolete: bool = False
    replaced_by: str | None = None
    inchikey: str | None = None
    cas: list[str] = field(default_factory=list)
    alt_ids: list[str] = field(default_factory=list)


class ChEBI:
    def __init__(self, terms: dict[str, Term]):
        self.terms = terms
        self.by_name: dict[str, str] = {}
        for t in terms.values():
            if not t.obsolete and t.name:
                self.by_name.setdefault(t.name.lower(), t.id)
        self._ancestors = lru_cache(maxsize=None)(self._ancestors_uncached)
        self.by_inchikey: dict[str, list[str]] = {}
        self.by_skeleton: dict[str, list[str]] = {}
        self.by_cas: dict[str, list[str]] = {}
        self.alt_to_primary: dict[str, str] = {}
        for t in terms.values():
            if t.obsolete:
                continue
            if t.inchikey:
                self.by_inchikey.setdefault(t.inchikey, []).append(t.id)
                self.by_skeleton.setdefault(t.inchikey[:14], []).append(t.id)
            for c in t.cas:
                self.by_cas.setdefault(c, []).append(t.id)
            for a in t.alt_ids:
                self.alt_to_primary[a] = t.id

    # ------------------------------------------------------------------ parsing
    @classmethod
    def from_obo(cls, path: Path) -> ChEBI:
        opener = gzip.open if str(path).endswith(".gz") else open
        terms: dict[str, Term] = {}
        cur: Term | None = None
        in_term = False
        with opener(path, "rt", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if line.startswith("["):
                    if cur is not None:
                        terms[cur.id] = cur
                    cur, in_term = None, line == "[Term]"
                    continue
                if not in_term or not line or ":" not in line:
                    continue
                key, _, val = line.partition(": ")
                val = val.split(" ! ")[0].strip()
                if key == "id":
                    cur = Term(id=val)
                elif cur is None:
                    continue
                elif key == "name":
                    cur.name = val
                elif key == "is_a":
                    cur.is_a.append(val.split()[0])
                elif key == "relationship":
                    parts = val.split()
                    if len(parts) >= 2 and parts[0] in ("has_role", "RO:0000087"):
                        cur.has_role.append(parts[1])
                elif key == "alt_id":
                    cur.alt_ids.append(val)
                elif key == "xref" and val.startswith("cas:"):
                    cur.cas.append(val[4:].split()[0])
                elif key == "property_value" and "inchi_key_string" in val:
                    parts = val.split('"')
                    if len(parts) >= 2:
                        cur.inchikey = parts[1]
                elif key == "is_obsolete" and val == "true":
                    cur.obsolete = True
                elif key == "replaced_by":
                    cur.replaced_by = val
        if cur is not None:
            terms[cur.id] = cur
        return cls(terms)

    # ------------------------------------------------------------------ queries
    def _ancestors_uncached(self, term_id: str) -> frozenset[str]:
        out: set[str] = set()
        stack = list(self.terms.get(term_id, Term(term_id)).is_a)
        while stack:
            t = stack.pop()
            if t in out:
                continue
            out.add(t)
            stack.extend(self.terms.get(t, Term(t)).is_a)
        return frozenset(out)

    def ancestors(self, term_id: str) -> frozenset[str]:
        return self._ancestors(term_id)

    def current_id(self, term_id: str) -> str:
        """Map alt IDs to their primary ID and follow ``replaced_by`` for obsolete IDs."""
        term_id = self.alt_to_primary.get(term_id, term_id)
        seen = set()
        while term_id in self.terms and self.terms[term_id].obsolete:
            nxt = self.terms[term_id].replaced_by
            if not nxt or nxt in seen:
                break
            seen.add(term_id)
            term_id = nxt
        return term_id

    def roles(self, term_id: str) -> set[str]:
        """Roles asserted directly on the term."""
        return set(self.terms.get(self.current_id(term_id), Term(term_id)).has_role)

    def role_id(self, name: str) -> str | None:
        return self.by_name.get(name.lower())

    def pesticide_classes(self, term_id: str) -> list[str]:
        """Named agrochemical classes (herbicide, fungicide, ...) asserted for the term.

        Returns ``["agrochemical"]`` when the term has an agrochemical or pesticide role
        that is not one of the named classes, and ``[]`` when it has none.
        """
        umbrella = {self.role_id(n) for n in ("pesticide", "agrochemical",
                                              "plant growth regulator")} - {None}
        named = {n: self.role_id(n) for n in PESTICIDE_CLASSES}
        out: set[str] = set()
        any_agro = False
        for r in self.roles(term_id):
            lineage = {r, *self.ancestors(r)}
            if lineage & umbrella:
                any_agro = True
                for n, nid in named.items():
                    if nid and nid in lineage:
                        out.add(n)
        if any_agro and not out:
            out.add("agrochemical")
        return sorted(out)
