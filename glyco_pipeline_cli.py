"""Command line interface for the glycan pipeline."""

import argparse
from pathlib import Path
from typing import List

from microberx.glyco_pipeline import (
    GlycanInputHandler,
    RuleGenerator,
    GlycanPredictor,
    EnzymeAnnotator,
    PathwayMapper,
    OutputGenerator,
    DataIntegrator,
)


def load_list(path: str) -> List[str]:
    try:
        with open(path) as fh:
            return [line.strip() for line in fh if line.strip()]
    except OSError:
        return []


def main() -> None:
    ap = argparse.ArgumentParser(description="Run glycan prediction pipeline")
    ap.add_argument("input", help="GlyTouCan ID, SMILES or file path")
    ap.add_argument("--format", dest="fmt", help="Input format")
    ap.add_argument("--metabolites", help="File of detected metabolites")
    ap.add_argument("--species", help="File of detected species")
    args = ap.parse_args()

    handler = GlycanInputHandler()
    glycan, mol = handler.load_input(args.input, fmt=args.fmt)

    rg = RuleGenerator()
    rules = rg.generate_rules([])

    predictor = GlycanPredictor(rules)
    graph = predictor.run(mol)

    annotator = EnzymeAnnotator()
    rule_species = {}
    for rule in rules:
        info = annotator.annotate_rule(rule.rule_id, rule.ec)
        rule_species[rule.rule_id] = info.get("species", [])

    mapper = PathwayMapper()

    metabolites = set(load_list(args.metabolites)) if args.metabolites else set()
    species = set(load_list(args.species)) if args.species else set()
    integrator = DataIntegrator(metabolites, species)

    product_scores = {k: v for k, v in predictor.seen.items()}
    adjusted_scores = integrator.adjust_scores(product_scores)

    outgen = OutputGenerator()
    report = outgen.report(graph)
    print(report)


if __name__ == "__main__":
    main()
