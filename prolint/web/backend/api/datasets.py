"""Dataset management API endpoints.

This module provides REST endpoints for uploading, listing, retrieving,
and deleting molecular dynamics datasets. Each dataset consists of a
topology file and optional trajectory file.
"""

from fastapi import APIRouter, HTTPException, UploadFile, File, Form
from typing import List, Optional

from prolint.web.backend.models.schemas import DatasetInfo
from prolint.web.backend.services import DatasetService

router = APIRouter()


@router.post("/upload", response_model=DatasetInfo)
async def upload_dataset(
    topology_file: UploadFile = File(...),
    trajectory_file: Optional[UploadFile] = File(None),
    name: Optional[str] = Form(None),
):
    """Upload a new molecular dynamics dataset.

    Parameters
    ----------
    topology_file : UploadFile
        Topology file (PDB, GRO, PSF, etc.).
    trajectory_file : UploadFile, optional
        Trajectory file (XTC, TRR, DCD, etc.).
    name : str, optional
        Custom name for the dataset.

    Returns
    -------
    DatasetInfo
        Metadata about the uploaded dataset.

    Raises
    ------
    HTTPException
        400 if file format is invalid, 500 for internal errors.
    """
    try:
        dataset_id = await DatasetService.load_from_upload(
            topology_file=topology_file.file,
            topology_filename=topology_file.filename or "topology",
            trajectory_file=trajectory_file.file if trajectory_file else None,
            trajectory_filename=trajectory_file.filename if trajectory_file else None,
            name=name,
        )

        metadata = DatasetService.get_dataset_info(dataset_id)
        return DatasetInfo(**metadata)

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Internal error: {str(e)}")


@router.get("/", response_model=List[DatasetInfo])
async def list_datasets():
    """List all available datasets.

    Returns
    -------
    list of DatasetInfo
        Metadata for all datasets in storage.

    Raises
    ------
    HTTPException
        500 for internal errors.
    """
    try:
        datasets = DatasetService.list_datasets()
        return [DatasetInfo(**ds) for ds in datasets]

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Internal error: {str(e)}")


@router.get("/{dataset_id}", response_model=DatasetInfo)
async def get_dataset(dataset_id: str):
    """Get metadata for a specific dataset.

    Parameters
    ----------
    dataset_id : str
        Unique identifier for the dataset.

    Returns
    -------
    DatasetInfo
        Metadata about the requested dataset.

    Raises
    ------
    HTTPException
        404 if dataset not found, 500 for internal errors.
    """
    try:
        metadata = DatasetService.get_dataset_info(dataset_id)

        if metadata is None:
            raise HTTPException(
                status_code=404, detail=f"Dataset not found: {dataset_id}"
            )

        return DatasetInfo(**metadata)

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Internal error: {str(e)}")


@router.delete("/{dataset_id}")
async def delete_dataset(dataset_id: str):
    """Delete a dataset from storage.

    Parameters
    ----------
    dataset_id : str
        Unique identifier for the dataset to delete.

    Returns
    -------
    dict
        Success message confirming deletion.

    Raises
    ------
    HTTPException
        404 if dataset not found, 500 for internal errors.
    """
    try:
        success = DatasetService.delete_dataset(dataset_id)

        if not success:
            raise HTTPException(
                status_code=404, detail=f"Dataset not found: {dataset_id}"
            )

        return {"message": f"Dataset {dataset_id} deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Internal error: {str(e)}")
