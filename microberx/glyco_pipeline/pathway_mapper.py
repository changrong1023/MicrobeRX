"""Map metabolites to known pathways."""

from __future__ import annotations

from typing import Dict, List


class PathwayMapper:
    def __init__(self, mapping: Dict[str, List[str]] | None = None) -> None:
        default = {
            "Glucose": ["Glycolysis", "Starch and sucrose metabolism"],
            "Galactose": ["Galactose metabolism"],
            "N-Acetylglucosamine": ["Amino sugar metabolism"],
        }
        self.mapping = mapping or default

    def map_compound(self, name: str) -> List[str]:
        return self.mapping.get(name, [])

    @classmethod
    def from_json(cls, path: str) -> "PathwayMapper":
        import json

        with open(path) as fh:
            data = json.load(fh)
        return cls(data)
