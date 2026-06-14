"""Dashboard API endpoints for ProLint web interface.

This module provides REST endpoints for computing and retrieving
biomolecular interaction analyses, including density maps, timeseries,
kinetics, and structure exports.
"""

from fastapi import APIRouter, HTTPException
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel
from typing import Optional
import numpy as np
import tempfile
import os

from prolint.web.backend.services.interaction_service import InteractionService
from prolint.web.backend.services.storage_service import storage

router = APIRouter()


def get_result_or_404(result_id: str):
    """Get result from storage or raise 404.

    Parameters
    ----------
    result_id : str
        Unique identifier for the computation result.

    Returns
    -------
    InteractionResult
        The stored computation result.

    Raises
    ------
    HTTPException
        404 if result not found.
    """
    result = storage.get_result(result_id)
    if result is None:
        raise HTTPException(status_code=404, detail=f"Result not found: {result_id}")
    return result


def require_param(value, name: str):
    """Validate that a required parameter is present.

    Parameters
    ----------
    value : any
        The parameter value to check.
    name : str
        Parameter name for error message.

    Returns
    -------
    any
        The original value if not None.

    Raises
    ------
    HTTPException
        400 if value is None.
    """
    if value is None:
        raise HTTPException(status_code=400, detail=f"{name} parameter is required")
    return value


class DashboardComputeRequest(BaseModel):
    """Request model for interaction computation.

    Attributes
    ----------
    dataset_id : str
        ID of the uploaded dataset to analyze.
    query_selection : str
        MDAnalysis selection string for query group.
    database_selection : str
        MDAnalysis selection string for database group.
    cutoff : float
        Contact distance cutoff in Angstroms.
    start : int
        Starting frame index.
    stop : int, optional
        Ending frame index (None for last frame).
    step : int
        Frame step interval.
    units : str
        Time units for output ("ns", "us", "frames").
    normalize_by : str
        Normalization method ("counts", "frames", "time").
    selected_replica : str, optional
        Selected replica ID (e.g., 'A', 'B'). Required when multiple
        replicas with repeated residue IDs are detected.
    replica_info : list, optional
        Information about detected replicas from replica analysis.
    """

    dataset_id: str
    query_selection: str
    database_selection: str
    cutoff: float = 7.0
    start: int = 0
    stop: Optional[int] = None
    step: int = 1
    units: str = "ns"
    normalize_by: str = "counts"
    selected_replica: Optional[str] = None
    replica_info: Optional[list] = None


@router.post("/compute")
async def compute_interactions(request: DashboardComputeRequest):
    """Compute biomolecular interactions for a dataset.

    Parameters
    ----------
    request : DashboardComputeRequest
        Computation parameters including selections and cutoff.

    Returns
    -------
    dict
        Result ID, computation time, and status.

    Raises
    ------
    HTTPException
        400 if parameters are invalid.
    """
    try:
        result_id, computation_time = InteractionService.compute_all_granularities(
            dataset_id=request.dataset_id,
            groupA_selection=request.query_selection,
            groupB_selection=request.database_selection,
            cutoff=request.cutoff,
            start=request.start,
            stop=request.stop,
            step=request.step,
            units=request.units,
            normalize_by=request.normalize_by,
            selected_replica=request.selected_replica,
            replica_info=request.replica_info,
        )
        return {
            "result_id": result_id,
            "computation_time": computation_time,
            "status": "completed",
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


class ReplicaAnalysisRequest(BaseModel):
    """Request model for replica analysis."""

    dataset_id: str
    query_selection: str


@router.post("/analyze-replicas")
async def analyze_replicas(request: ReplicaAnalysisRequest):
    """Analyze replica structure in a query selection.

    Detects multiple replicas (e.g., protein replicates) and checks
    for repeated residue IDs to determine how they should be handled.

    Parameters
    ----------
    request : ReplicaAnalysisRequest
        Dataset ID and query selection string.

    Returns
    -------
    dict
        Replica analysis results including:
        - n_replicas: number of replicas detected
        - detection_method: how replicas were detected
        - has_repeated_resids: whether replicas have repeated residue IDs
        - replica_info: details about each replica

    Raises
    ------
    HTTPException
        400 if selection is invalid, 404 if dataset not found.
    """
    from prolint.web.backend.services import DatasetService

    try:
        result = DatasetService.analyze_query_replicas(
            dataset_id=request.dataset_id,
            query_selection=request.query_selection,
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/{result_id}/interactions")
async def get_interaction_data(result_id: str):
    """Get interaction summary data.

    Parameters
    ----------
    result_id : str
        Computation result ID.

    Returns
    -------
    dict
        Result ID and interaction summary data.

    Raises
    ------
    HTTPException
        404 if result not found.
    """
    result = get_result_or_404(result_id)
    return {"result_id": result_id, "interactions": result.get_interaction_data()}


@router.get("/{result_id}/density-map")
async def get_density_map(
    result_id: str,
    frame_start: int = 0,
    frame_end: int = None,
    bins: int = 50,
    database_types: str = None,
):
    """Get 2D density map data.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    frame_start : int, default=0
        Starting frame for density calculation.
    frame_end : int, optional
        Ending frame for density calculation.
    bins : int, default=50
        Number of bins for the density grid.
    database_types : str, optional
        Comma-separated database residue types to include.

    Returns
    -------
    dict
        Density map data with grid and values.

    Raises
    ------
    HTTPException
        404 if result not found.
    """
    result = get_result_or_404(result_id)
    selected_types = (
        [t.strip() for t in database_types.split(",") if t.strip()]
        if database_types
        else None
    )
    analysis = result.contacts.analyze(
        "density_map",
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=result.frame_step,
        bins=bins,
        database_types=selected_types,
    )
    return analysis.data


@router.get("/{result_id}/structure")
async def get_structure_with_metrics(
    result_id: str,
    metric: str = "occupancy",
    frame_idx: int = 0,
    database_type: str = None,
):
    """Get PDB structure with metric values in B-factor column.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    metric : str, default="occupancy"
        Metric to encode in B-factors ("occupancy", "mean", "max").
    frame_idx : int, default=0
        Frame index for atomic coordinates.
    database_type : str
        Database residue type for metric calculation.

    Returns
    -------
    PlainTextResponse
        PDB file content with metrics in B-factor column.

    Raises
    ------
    HTTPException
        400 if database_type missing, 404 if result/atoms not found.
    """
    result = get_result_or_404(result_id)
    require_param(database_type, "database_type")

    universe = result.universe
    universe.trajectory[frame_idx]
    query_atoms = universe.query.atoms

    if len(query_atoms) == 0:
        raise HTTPException(status_code=404, detail="No query atoms found")

    # Compute metrics
    metric_results = result.contacts.compute_metric(
        metric, target_resname=database_type
    )
    residue_metrics = {
        int(resid): data[database_type]["global"]
        for resid, data in metric_results.items()
        if database_type in data
    }

    # Normalize to 0-100 for B-factors
    if residue_metrics:
        max_val = max(residue_metrics.values()) or 1
        residue_metrics = {k: (v / max_val) * 100 for k, v in residue_metrics.items()}

    # Set B-factors
    if not hasattr(universe.atoms, "tempfactors"):
        from MDAnalysis.core.topologyattrs import Tempfactors

        universe.add_TopologyAttr(Tempfactors(np.zeros(len(universe.atoms))))

    for atom in query_atoms:
        atom.tempfactor = residue_metrics.get(int(atom.residue.resid), 0.0)

    # Write PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        tmp_filename = tmp.name

    try:
        query_atoms.write(tmp_filename)
        with open(tmp_filename, "r") as f:
            pdb_content = f.read()
    finally:
        if os.path.exists(tmp_filename):
            os.remove(tmp_filename)

    return PlainTextResponse(content=pdb_content, media_type="chemical/x-pdb")


@router.get("/{result_id}/logoplot")
async def get_logoplot_data(
    result_id: str, metric: str = "occupancy", database_type: str = None
):
    """Get per-residue metric data for logo plot.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    metric : str, default="occupancy"
        Metric to calculate per residue.
    database_type : str
        Database residue type for metric calculation.

    Returns
    -------
    dict
        Residue data with resid, resname, chainID, and metric value.

    Raises
    ------
    HTTPException
        400 if database_type missing, 404 if result not found.
    """
    result = get_result_or_404(result_id)
    require_param(database_type, "database_type")

    metric_results = result.contacts.compute_metric(
        metric, target_resname=database_type
    )

    # Build residue data from the query (already filtered to selected replica)
    residue_data = []
    for res in result.universe.query.residues:
        chain_id = getattr(res, "segid", "A") if hasattr(res, "segid") else "A"
        residue_data.append({
            "resid": int(res.resid),
            "resname": res.resname,
            "chainID": result.selected_replica or chain_id,
            "value": (
                metric_results.get(int(res.resid), {})
                .get(database_type, {})
                .get("global", 0.0)
            ),
        })
    residue_data.sort(key=lambda x: x["resid"])

    return {"residues": residue_data}


@router.get("/{result_id}/shared-contacts")
async def get_shared_contacts_data(result_id: str, database_type: str = None):
    """Get shared contacts correlation matrix.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    database_type : str
        Database residue type for analysis.

    Returns
    -------
    dict
        Labels (residue IDs) and correlation matrix.

    Raises
    ------
    HTTPException
        400 if database_type missing, 404 if result not found.
    """
    result = get_result_or_404(result_id)
    require_param(database_type, "database_type")
    analysis = result.contacts.analyze("shared_contacts", database_type=database_type)
    return {"labels": analysis.data["labels"], "matrix": analysis.data["matrix"]}


@router.get("/{result_id}/timeseries")
async def get_timeseries_data(
    result_id: str, database_type: str = None, query_residues: str = None
):
    """Get contact timeseries for selected residues.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    database_type : str
        Database residue type for contact counting.
    query_residues : str
        Comma-separated query residue IDs.

    Returns
    -------
    dict
        Query residues, frames, and contact counts per residue.

    Raises
    ------
    HTTPException
        400 if parameters missing/invalid, 404 if result not found.
    """
    result = get_result_or_404(result_id)
    require_param(database_type, "database_type")
    require_param(query_residues, "query_residues")

    resid_list = [int(r.strip()) for r in query_residues.split(",") if r.strip()]
    if not resid_list:
        raise HTTPException(status_code=400, detail="No valid query residues provided")

    analysis = result.contacts.analyze(
        "timeseries",
        database_type=database_type,
        query_residues=resid_list,
        frame_start=result.frame_start,
        frame_end=result.frame_end,
        frame_step=result.frame_step,
    )
    return {
        "query_residues": analysis.data["query_residues"],
        "frames": analysis.data["frames"],
        "contact_counts": analysis.data["contact_counts"],
    }


@router.get("/{result_id}/distance-time")
async def get_distance_time(
    result_id: str, query_residue: int = None, database_residue: int = None
):
    """Get distance over time between a residue pair.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    query_residue : int
        Query residue ID.
    database_residue : int
        Database residue ID.

    Returns
    -------
    dict
        Frames, distances, min distances, contact frames, and positions.

    Raises
    ------
    HTTPException
        400 if parameters missing, 404 if result/residues not found.
    """
    result = get_result_or_404(result_id)
    require_param(query_residue, "query_residue")
    require_param(database_residue, "database_residue")

    try:
        analysis = result.contacts.analyze(
            "distances",
            query_residue=query_residue,
            database_residue=database_residue,
            frame_start=result.frame_start,
            frame_end=result.frame_end,
            frame_step=result.frame_step,
        )
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))

    return {
        "frames": analysis.data["frames"],
        "distances": analysis.data["distances"],
        "min_distances": analysis.data.get("min_distances", []),
        "contact_frames": analysis.data["contact_frames"],
        "positions": analysis.data.get("positions", {}),
    }


@router.get("/{result_id}/structure-interaction")
async def get_structure_interaction(
    result_id: str,
    query_residue: int = None,
    database_residue: int = None,
    frame_idx: int = 0,
):
    """Get PDB structure for interaction visualization.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    query_residue : int
        Query residue ID to include.
    database_residue : int, optional
        Database residue ID to include.
    frame_idx : int, default=0
        Frame index for atomic coordinates.

    Returns
    -------
    PlainTextResponse
        PDB file content for visualization.

    Raises
    ------
    HTTPException
        400 if query_residue missing, 404 if result/residues not found.
    """
    result = get_result_or_404(result_id)
    require_param(query_residue, "query_residue")

    universe = result.universe
    frame_idx = max(0, min(frame_idx, universe.trajectory.n_frames - 1))
    universe.trajectory[frame_idx]

    query_res_atoms = universe.query.select_atoms(f"resid {query_residue}")
    if len(query_res_atoms) == 0:
        raise HTTPException(
            status_code=404, detail=f"Query residue {query_residue} not found"
        )

    query_atoms = universe.query.atoms
    if database_residue is not None:
        db_atoms = universe.select_atoms(f"resid {database_residue}")
        if len(db_atoms) == 0:
            raise HTTPException(
                status_code=404, detail=f"Database residue {database_residue} not found"
            )
        combined = query_atoms | db_atoms
    else:
        combined = query_atoms

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        tmp_filename = tmp.name

    try:
        combined.write(tmp_filename)
        with open(tmp_filename, "r") as f:
            pdb_content = f.read()
    finally:
        if os.path.exists(tmp_filename):
            os.remove(tmp_filename)

    return PlainTextResponse(content=pdb_content, media_type="chemical/x-pdb")


@router.get("/{result_id}/residue-timeseries")
async def get_residue_timeseries(
    result_id: str, query_residue: int = None, database_type: str = None
):
    """Get per-database-residue contact timeline for a query residue.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    query_residue : int
        Query residue ID to analyze.
    database_type : str
        Database residue type to track.

    Returns
    -------
    dict
        Database IDs, frames, and contact matrix.

    Raises
    ------
    HTTPException
        400 if parameters missing, 404 if result not found.
    """
    result = get_result_or_404(result_id)
    require_param(query_residue, "query_residue")
    require_param(database_type, "database_type")

    analysis = result.contacts.analyze(
        "database_contacts",
        query_residue=query_residue,
        database_type=database_type,
        frame_start=result.frame_start,
        frame_end=result.frame_end,
        frame_step=result.frame_step,
        top_n=50,
    )

    return {
        "database_ids": analysis.data["database_ids"],
        "frames": analysis.data["frames"],
        "contact_matrix": analysis.data["contact_matrix"],
        "total_database_ids": analysis.data["total_database_ids"],
    }


@router.get("/{result_id}/atom-distances")
async def get_atom_distances(
    result_id: str,
    query_residue: int = None,
    database_residue: int = None,
    frame_idx: int = None,
):
    """Get atom-atom distance matrix between two residues.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    query_residue : int
        Query residue ID.
    database_residue : int
        Database residue ID.
    frame_idx : int
        Frame index for distance calculation.

    Returns
    -------
    dict
        Frame, atom labels, distance matrix, and min/max distances.

    Raises
    ------
    HTTPException
        400 if parameters missing, 404 if result/residues not found.
    """
    result = get_result_or_404(result_id)
    require_param(query_residue, "query_residue")
    require_param(database_residue, "database_residue")
    require_param(frame_idx, "frame_idx")

    try:
        analysis = result.contacts.analyze(
            "atom_distances",
            query_residue=query_residue,
            database_residue=database_residue,
            frame_idx=frame_idx,
        )
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))

    return {
        "frame": analysis.data["frame"],
        "query_atoms": analysis.data["query_atoms"],
        "database_atoms": analysis.data["database_atoms"],
        "distance_matrix": analysis.data["distance_matrix"],
        "min_distance": analysis.data["min_distance"],
        "max_distance": analysis.data["max_distance"],
    }


@router.get("/{result_id}/kinetics")
async def get_kinetics_data(
    result_id: str,
    query_residue: int = None,
    database_residue: int = None,
    database_type: str = None,
    mode: str = "individual",
):
    """Get kinetics analysis data including survival curves.

    Parameters
    ----------
    result_id : str
        Computation result ID.
    query_residue : int
        Query residue ID.
    database_residue : int, optional
        Database residue ID (required for individual mode).
    database_type : str, optional
        Database residue type (required for accumulated mode).
    mode : str, default="individual"
        Analysis mode ("individual" or "accumulated").

    Returns
    -------
    dict
        Kinetics data with survival curves and fit parameters.

    Raises
    ------
    HTTPException
        400 if required parameters missing for mode, 404 if not found.
    """
    result = get_result_or_404(result_id)
    require_param(query_residue, "query_residue")

    if mode == "individual" and database_residue is None:
        raise HTTPException(
            status_code=400, detail="database_residue required for individual mode"
        )
    if mode == "accumulated" and database_type is None:
        raise HTTPException(
            status_code=400, detail="database_type required for accumulated mode"
        )

    try:
        analysis = result.contacts.analyze(
            "kinetics",
            query_residue=query_residue,
            database_residue=database_residue,
            database_type=database_type,
            mode=mode,
            fit_survival=True,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    return analysis.data
