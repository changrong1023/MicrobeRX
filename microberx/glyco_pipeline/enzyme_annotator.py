"""Annotate reaction rules with enzyme and species information."""

from __future__ import annotations

import json
from typing import Dict, List, Optional


class EnzymeAnnotator:
    def __init__(self, db: Dict[str, Dict[str, List[str]]] | None = None) -> None:
        """db maps rule_id or EC to species list"""
        default = {
            "3.2.1.x": {"species": ["Escherichia coli", "Bifidobacterium longum"]}
        }
        self.db = db or default

    def annotate_rule(self, rule_id: str, ec: Optional[str] = None) -> Dict[str, List[str]]:
        """Return species list for a rule."""
        info = self.db.get(rule_id) or self.db.get(ec or "") or {}
        return info

    @classmethod
    def from_json(cls, path: str) -> "EnzymeAnnotator":
        with open(path) as fh:
            data = json.load(fh)
        return cls(data)
