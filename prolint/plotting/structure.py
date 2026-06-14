"""Structure export module.

This module provides functions for exporting contact metrics
to PDB files for visualization in molecular viewers.
"""

from typing import Optional, Literal
import tempfile


def write_pdb(
    contacts,
    metric: Literal["mean", "max", "sum", "occupancy"] = "occupancy",
    target_resname: Optional[str] = None,
    filename: Optional[str] = None,
    frame: int = 0,
) -> str:
    """Write contact metrics to a PDB file for visualization.

    Exports query atoms to a PDB file with metric values stored in the
    B-factor column for coloring in molecular viewers.

    Parameters
    ----------
    contacts : ComputedContacts
        Computed contact data.
    metric : {"mean", "max", "sum", "occupancy"}, default="occupancy"
        Metric to write to B-factor column.
    target_resname : str, optional
        Filter by database residue name (e.g., "CHOL").
    filename : str, optional
        Output filename. If None, creates a temporary file.
    frame : int, default=0
        Trajectory frame to use for coordinates.

    Returns
    -------
    str
        Path to the written PDB file.

    Examples
    --------
    >>> from prolint.plotting import write_pdb
    >>> pdb_path = write_pdb(contacts, metric="occupancy")
    >>> # Open in PyMOL/VMD and color by B-factor
    """
    universe = contacts.provider.query.universe
    query = contacts.provider.query

    # Compute metric values
    metrics = contacts.compute_metric(metric, target_resname=target_resname)

    # Build resid -> value mapping (use global value for each residue)
    bfactors = {}
    for resid, db_data in metrics.items():
        # Sum across all database types if target_resname is None
        if target_resname:
            bfactors[resid] = db_data.get(target_resname, {}).get("global", 0.0)
        else:
            # Average across all database types
            values = [d["global"] for d in db_data.values()]
            bfactors[resid] = sum(values) / len(values) if values else 0.0

    # Set frame
    universe.trajectory[frame]

    # Generate filename if not provided
    if filename is None:
        fd, filename = tempfile.mkstemp(suffix=".pdb", prefix="prolint_")
        import os

        os.close(fd)

    # Write PDB
    atoms = query.atoms
    with open(filename, "w") as f:
        f.write(f"REMARK    Metric: {metric}\n")
        if target_resname:
            f.write(f"REMARK    Database: {target_resname}\n")
        f.write(f"REMARK    Frame: {frame}\n")

        if universe.dimensions is not None:
            dims = universe.dimensions
            f.write(
                f"CRYST1{dims[0]:9.3f}{dims[1]:9.3f}{dims[2]:9.3f}"
                f"{dims[3]:7.2f}{dims[4]:7.2f}{dims[5]:7.2f} P 1           1\n"
            )

        for i, atom in enumerate(atoms):
            serial = (i + 1) % 100000
            name = atom.name[:4].ljust(4)
            resname = atom.resname[:3].ljust(3)
            chain = getattr(atom, "chainID", " ") or " "
            resid = atom.resid % 10000
            x, y, z = atom.position
            occupancy = getattr(atom, "occupancy", 1.0)
            tempfactor = bfactors.get(atom.resid, 0.0)
            element = getattr(atom, "element", atom.name[0])[:2].rjust(2)
            f.write(
                f"ATOM  {serial:5d} {name}{' '}{resname} {chain}{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempfactor:6.2f}"
                f"          {element}\n"
            )

        f.write("END\n")

    return filename
