"""Multi-step glycan transformation predictor."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, List, Tuple

from rdkit import Chem

from .reaction_rules import ReactionRule

logger = logging.getLogger(__name__)


class GlycanPredictor:
    """Iteratively apply reaction rules to a glycan molecule."""

    def __init__(self, rules: List[ReactionRule], max_steps: int = 5, score_threshold: float = 0.1) -> None:
        self.rules = rules
        self.max_steps = max_steps
        self.score_threshold = score_threshold
        self.graph: Dict[str, List[Tuple[str, List[str], float]]] = defaultdict(list)
        self.seen: Dict[str, float] = {}

    def _mol_key(self, mol: Chem.Mol) -> str:
        return Chem.MolToSmiles(mol, canonical=True)

    def run(self, start_mol: Chem.Mol) -> Dict[str, List[Tuple[str, List[str]]]]:
        current = [(start_mol, 1.0)]
        self.seen[self._mol_key(start_mol)] = 1.0
        step = 0
        while current and step < self.max_steps:
            next_gen = []
            for mol, score in current:
                mkey = self._mol_key(mol)
                for rule in self.rules:
                    new_score = score * rule.score
                    if new_score < self.score_threshold:
                        continue
                    for prods in rule.apply(mol):
                        prod_keys = []
                        for p in prods:
                            key = self._mol_key(p)
                            prod_keys.append(key)
                            prev_score = self.seen.get(key, 0.0)
                            if new_score > prev_score:
                                self.seen[key] = new_score
                                next_gen.append((p, new_score))
                        self.graph[mkey].append((rule.rule_id, prod_keys, new_score))
            current = next_gen
            step += 1
        return self.graph
