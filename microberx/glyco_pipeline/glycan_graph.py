"""Simple graph representation for glycans."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass
class MonosaccharideNode:
    """A monosaccharide residue in a glycan."""

    id: str
    name: str
    anomeric: str | None = None
    ring_size: str | None = None
    mods: Dict[str, str] | None = None
    children: List[Tuple[str, int, int]] = field(default_factory=list)


@dataclass
class GlycanGraph:
    """Graph model consisting of monosaccharide nodes."""

    nodes: Dict[str, MonosaccharideNode] = field(default_factory=dict)
    root: str | None = None

    def add_node(self, node: MonosaccharideNode) -> None:
        self.nodes[node.id] = node

    def add_edge(self, parent_id: str, child_id: str, parent_pos: int, child_pos: int) -> None:
        parent = self.nodes[parent_id]
        parent.children.append((child_id, parent_pos, child_pos))
