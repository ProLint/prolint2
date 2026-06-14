"""
Pydantic Schemas for API Requests and Responses
"""

from pydantic import BaseModel, Field
from typing import Optional


class DatasetInfo(BaseModel):
    """Dataset information."""

    id: str = Field(..., description="Unique dataset ID")
    name: str = Field(..., description="Dataset name")
    topology_file: Optional[str] = Field(None, description="Path to topology file")
    trajectory_file: Optional[str] = Field(None, description="Path to trajectory file")
    n_frames: Optional[int] = Field(None, description="Total frames in trajectory")
    n_atoms: Optional[int] = Field(None, description="Total atoms")
    n_residues: Optional[int] = Field(None, description="Total residues")
    status: str = Field(default="ready", description="Dataset status")
