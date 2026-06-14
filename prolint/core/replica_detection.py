"""Replica detection utilities for ProLint."""

from dataclasses import dataclass, field
from typing import Optional, List, Any
import logging

import numpy as np

logger = logging.getLogger(__name__)

SEGMENT_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


@dataclass
class ReplicaInfo:
    """Information about a single molecular replica."""
    replica_id: str
    n_residues: int
    first_resid: int
    last_resid: int


@dataclass
class ReplicaDetectionResult:
    """Result of replica detection analysis."""
    n_replicas: int
    detection_method: Optional[str]
    has_repeated_resids: bool
    replica_info: List[ReplicaInfo] = field(default_factory=list)
    fragments: List[Any] = field(default_factory=list)


def detect_replicas(query_atoms) -> ReplicaDetectionResult:
    """Detect replicas in a query atom selection.

    Tries bond-based connectivity first, then falls back to residue ID
    replication detection for systems without bond info.
    """
    atoms = query_atoms.atoms if hasattr(query_atoms, "atoms") else query_atoms

    if len(atoms) == 0:
        return ReplicaDetectionResult(0, None, False)

    fragments = None
    detection_method = None
    has_repeated_resids = False

    # Try bond-based fragment detection
    try:
        if len(atoms.bonds) > 0:
            fragments = list(atoms.fragments)
            if len(fragments) > 1:
                detection_method = "bond_connectivity"
                logger.info(f"Detected {len(fragments)} replicas by bond connectivity")
    except Exception:
        pass

    # Fallback: detect by replicated residue IDs
    if not fragments or len(fragments) <= 1:
        resid_fragments = _detect_fragments_by_resid(atoms)
        if len(resid_fragments) > 1:
            fragments = resid_fragments
            detection_method = "resid_replication"
            has_repeated_resids = True
            logger.info(f"Detected {len(fragments)} replicas by residue ID replication")

    # Single replica or no fragments
    if not fragments or len(fragments) <= 1:
        return ReplicaDetectionResult(1, None, False, [], [atoms] if len(atoms) > 0 else [])

    # Build replica info
    replica_info = [
        ReplicaInfo(
            replica_id=SEGMENT_IDS[i % len(SEGMENT_IDS)],
            n_residues=len(frag.residues),
            first_resid=int(frag.residues[0].resid) if len(frag.residues) > 0 else 0,
            last_resid=int(frag.residues[-1].resid) if len(frag.residues) > 0 else 0,
        )
        for i, frag in enumerate(fragments)
    ]

    return ReplicaDetectionResult(
        len(fragments), detection_method, has_repeated_resids, replica_info, fragments
    )


def _detect_fragments_by_resid(query_atoms) -> List:
    """Detect replicas by finding replicated residue ID sequences."""
    residues = query_atoms.residues
    if len(residues) == 0:
        return []

    resids = residues.resids
    unique_resids, counts = np.unique(resids, return_counts=True)

    if np.all(counts == 1):
        return []

    n_replicates = int(np.max(counts))

    # Build mapping of resid to residue indices
    resid_to_indices = {rid: [] for rid in unique_resids}
    for idx, rid in enumerate(resids):
        resid_to_indices[rid].append(idx)

    # Distribute residues across replicas
    fragment_atom_indices = [[] for _ in range(n_replicates)]
    for rid in unique_resids:
        for frag_idx, res_idx in enumerate(resid_to_indices[rid]):
            if frag_idx < n_replicates:
                fragment_atom_indices[frag_idx].extend(residues[res_idx].atoms.ix)

    # Convert to AtomGroups
    universe = query_atoms.universe
    return [universe.atoms[idx] for idx in fragment_atom_indices if idx]


def get_replica_atoms(replica_result: ReplicaDetectionResult, replica_id: str):
    """Get atoms for a specific replica by ID. Raises ValueError if not found."""
    for info, fragment in zip(replica_result.replica_info, replica_result.fragments):
        if info.replica_id == replica_id:
            return fragment

    available = [info.replica_id for info in replica_result.replica_info]
    raise ValueError(f"Replica '{replica_id}' not found. Available: {available}")
