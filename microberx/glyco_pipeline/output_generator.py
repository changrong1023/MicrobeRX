"""Generate textual report for predictions."""

from __future__ import annotations

from typing import Dict, List


class OutputGenerator:
    def __init__(self) -> None:
        pass

    def report(self, graph: Dict[str, List[tuple[str, List[str], float]]]) -> str:
        lines = ["# Prediction Results\n"]
        dot = ["digraph G {"]
        for parent, edges in graph.items():
            for rule_id, products, score in edges:
                prod_str = ", ".join(products)
                lines.append(f"{parent} --{rule_id} ({score:.2f})--> {prod_str}")
                for p in products:
                    dot.append(f'"{parent}" -> "{p}" [label="{rule_id} ({score:.2f})"];')
        dot.append("}")
        lines.append("\nGraphviz DOT:\n" + "\n".join(dot))
        return "\n".join(lines)
