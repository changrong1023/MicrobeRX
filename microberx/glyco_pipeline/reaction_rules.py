"""Reaction rule objects and generator for glycans."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterable, List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass
class ReactionRule:
    """Representation of a reaction rule."""

    rule_id: str
    reaction_smarts: str
    name: str | None = None
    ec: str | None = None
    enzyme: str | None = None
    score: float = 1.0

    def __post_init__(self) -> None:
        self.compiled = None
        try:
            self.compiled = AllChem.ReactionFromSmarts(self.reaction_smarts)
        except Exception:
            self.compiled = None

    def apply(self, mol: Chem.Mol) -> List[tuple[Chem.Mol, ...]]:
        if self.compiled is None:
            return []
        try:
            return list(self.compiled.RunReactants((mol,)))
        except Exception:
            return []


class RuleGenerator:
    """Generate reaction rules from data sources."""

    def __init__(self) -> None:
        self.rules: list[ReactionRule] = []

    def load_from_file(self, path: str) -> None:
        with open(path) as fh:
            data = json.load(fh)
        self.rules = [ReactionRule(**d) for d in data]

    def save_to_file(self, path: str) -> None:
        with open(path, "w") as fh:
            json.dump([r.__dict__ for r in self.rules], fh, indent=2)

    def generate_rules(self, sources: Iterable[str]) -> list[ReactionRule]:
        """Generate a minimal set of glycosidic bond hydrolysis rules."""
        self.rules = []
        idx = 1
        for conf in ["alpha", "beta"]:
            for pos in range(2, 7):
                smarts = (
                    f"[C@H:1]-O[C@@H:2]>>[C@H:1]-O.[C@@H:2]-O"  # simplified
                )
                rule = ReactionRule(
                    rule_id=f"R{idx:03d}",
                    reaction_smarts=smarts,
                    name=f"{conf}-1,{pos} cleavage",
                    ec="3.2.1.x",
                    enzyme="glycosidase",
                    score=1.0,
                )
                self.rules.append(rule)
                idx += 1
        return self.rules
