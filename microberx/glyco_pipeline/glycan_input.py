"""Glycan input handler for GlycoMicrobeRX pipeline."""

from __future__ import annotations

import logging
from typing import Optional, Tuple

try:
    from glypy import Glycan
except ImportError:  # pragma: no cover - optional dependency
    Glycan = None

from rdkit import Chem
import requests

from .glycan_graph import GlycanGraph, MonosaccharideNode


logger = logging.getLogger(__name__)


class GlycanInputHandler:
    """Load and parse glycan inputs from IDs or files."""

    def __init__(self, glytoucan_endpoint: str | None = None):
        self.glytoucan_endpoint = (
            glytoucan_endpoint
            or "https://sparqlist.glyconavi.org/api/GTC-WURCS?id={}"
        )

    def load_input(
        self, input_str: str, fmt: Optional[str] = None
    ) -> Tuple[object, Chem.Mol | None]:
        """Load a glycan structure.

        Parameters
        ----------
        input_str:
            GlyTouCan ID, file path or direct structure string.
        fmt:
            Explicit format: ``"wurcs"``, ``"smiles"`` or ``"glycoct"``.

        Returns
        -------
        glycan_obj:
            Parsed glycan object (glypy ``Glycan`` when available).
        rdkit_mol:
            Equivalent RDKit molecule if possible.
        """
        if fmt is None:
            fmt = self._guess_format(input_str)
        logger.info("Detected input format: %s", fmt)

        if fmt == "glytoucan":
            wurcs = self._fetch_wurcs(input_str)
            fmt = "wurcs"
            input_str = wurcs

        if fmt == "wurcs" and Glycan is not None:
            glycan = Glycan.from_wurcs(input_str)
            graph = self._build_graph(glycan)
            smiles = glycan.to_smiles()
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            return graph, mol

        if fmt == "smiles":
            mol = Chem.MolFromSmiles(input_str)
            if mol is None:
                raise ValueError("SMILES parsing failed")
            graph = GlycanGraph()
            node = MonosaccharideNode(id="0", name="SMILES")
            graph.add_node(node)
            graph.root = "0"
            return graph, mol

        if fmt == "glycoct" and Glycan is not None:
            glycan = Glycan.from_glycoct(input_str)
            graph = self._build_graph(glycan)
            smiles = glycan.to_smiles()
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            return graph, mol

        raise ValueError(f"Unsupported or unknown format: {fmt}")

    def _guess_format(self, s: str) -> str:
        if s.startswith("WURCS="):
            return "wurcs"
        if s.startswith("G") and len(s) == 8:
            return "glytoucan"
        if all(c.isalnum() for c in s) and len(s) == 8:
            return "glytoucan"
        if "." in s or "/" in s:
            # maybe filepath
            try:
                with open(s) as fh:
                    text = fh.read().strip()
                return self._guess_format(text)
            except OSError:
                pass
        return "smiles"

    def _fetch_wurcs(self, glytoucan_id: str) -> str:
        url = self.glytoucan_endpoint.format(glytoucan_id)
        logger.info("Fetching %s from GlyTouCan", glytoucan_id)
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        return data.get("wurcs", "")

    def _build_graph(self, glycan: Glycan) -> GlycanGraph:
        """Convert a ``glypy`` Glycan into ``GlycanGraph``."""
        graph = GlycanGraph()
        node_map: dict[int, str] = {}
        for res in glycan.iternodes():
            node_id = str(res.id)
            node = MonosaccharideNode(
                id=node_id,
                name=getattr(res, "name", "?"),
                anomeric=getattr(res, "anomeric_state", None),
                ring_size=getattr(res, "ring_type", None),
            )
            graph.add_node(node)
            node_map[res.id] = node_id
        if glycan.root:
            graph.root = node_map[glycan.root.id]
        for res in glycan.iternodes():
            if res.parent is not None:
                parent_id = node_map[res.parent.id]
                child_id = node_map[res.id]
                parent_pos = getattr(res, "parent_position", 1) or 1
                child_pos = getattr(res, "substituent_position", 1) or 1
                graph.add_edge(parent_id, child_id, parent_pos, child_pos)
        return graph
