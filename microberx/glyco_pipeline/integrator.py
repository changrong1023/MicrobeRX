"""Integrate metabolomic and genomic data."""

from __future__ import annotations

from typing import Dict, Iterable, Set


class DataIntegrator:
    def __init__(self, metabolites: Set[str] | None = None, species: Set[str] | None = None) -> None:
        self.metabolites = metabolites or set()
        self.species = species or set()

    def filter_metabolites(self, predicted: Iterable[str]) -> Dict[str, bool]:
        return {m: m in self.metabolites for m in predicted}

    def filter_species(self, rule_to_species: Dict[str, Iterable[str]]) -> Dict[str, bool]:
        out = {}
        for rule, sp in rule_to_species.items():
            out[rule] = bool(self.species.intersection(set(sp)))
        return out

    def adjust_scores(self, product_scores: Dict[str, float]) -> Dict[str, float]:
        """Boost scores if detected in datasets, penalize otherwise."""
        adjusted = {}
        for prod, score in product_scores.items():
            factor = 1.0
            if self.metabolites:
                if prod in self.metabolites:
                    factor *= 1.2
                else:
                    factor *= 0.8
            adjusted[prod] = score * factor
        return adjusted
