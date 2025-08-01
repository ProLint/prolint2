r"""UFCC sample data
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2025
:Copyright: MIT License
"""

import os
import time

import requests
import tqdm


def file_download(url: str, filename: str, desc: str, max_retries: int = 3):
    """Download a file from URL with progress bar and error handling.
    
    Args:
        url: URL to download from
        filename: Local filename to save to
        desc: Description for progress bar
        max_retries: Maximum number of retry attempts
    """
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    for attempt in range(max_retries):
        try:
            with open(filename, "wb") as f:
                with requests.get(url, stream=True, timeout=30) as r:
                    r.raise_for_status()
                    total = int(r.headers.get("content-length", 0))

                    tqdm_params = {
                        "total": total,
                        "miniters": 1,
                        "unit": "B",
                        "unit_scale": True,
                        "unit_divisor": 1024,
                        "desc": desc,
                    }
                    with tqdm.tqdm(**tqdm_params) as pb:
                        for chunk in r.iter_content(chunk_size=8192):
                            pb.update(len(chunk))
                            f.write(chunk)
            return  # Success - exit function
        except (requests.RequestException, IOError) as e:
            if attempt == max_retries - 1:
                raise RuntimeError(f"Failed to download {url} after {max_retries} attempts: {e}")
            time.sleep(2 ** attempt)  # Exponential backoff


class _DataSampleBase:
    """Base class for sample data."""
    
    def __init__(self, dataset_name: str, coordinates_url: str, trajectory_url: str):
        self.dataset_name = dataset_name
        self.path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), f"data/{dataset_name}/"
        )
        self._ensure_files_exist(coordinates_url, trajectory_url)
        self._set_file_paths()
    
    def _ensure_files_exist(self, coordinates_url: str, trajectory_url: str):
        """Download files if they don't exist."""
        coord_file = os.path.join(self.path, "coordinates.gro")
        traj_file = os.path.join(self.path, "trajectory.xtc")
        
        if not os.path.isfile(coord_file) or not os.path.isfile(traj_file):
            if not os.path.isfile(coord_file):
                file_download(
                    coordinates_url,
                    coord_file,
                    f"Downloading {self.dataset_name} coordinates...",
                )
            if not os.path.isfile(traj_file):
                file_download(
                    trajectory_url,
                    traj_file,
                    f"Downloading {self.dataset_name} trajectory...",
                )
    
    def _set_file_paths(self):
        """Set file path attributes."""
        self.coordinates = os.path.join(self.path, "coordinates.gro")
        self.trajectory = os.path.join(self.path, "trajectory.xtc")
        self.contacts = os.path.join(self.path, "contacts.csv")


class GIRKDataSample(_DataSampleBase):
    """
    Class to add sample data in order to test the installation of the package

    It uses the coordinates and trajectory of a simulation of the GIRK channel using the MARTINI force field.

    Example
    -------
    Import the prolint2 library and use the sample data as follows::

        from prolint2 import PL2
        from prolint2.sampledata import GIRK

        target_system = PL2(GIRK.coordinates, GIRK.trajectory)

    """

    def __init__(self):
        super().__init__(
            dataset_name="GIRK",
            coordinates_url="https://prolint.github.io/ProLintData/GIRK/mid/coordinates.gro",
            trajectory_url="https://prolint.github.io/ProLintData/GIRK/mid/trajectory.xtc"
        )


class COX1DataSample(_DataSampleBase):
    """
    Class to add sample data in order to test the installation of the package

    It uses the coordinates and trajectory of a simulation of the COX1 enzyme using the MARTINI force field.

    Example
    -------
    Import the prolint2 library and use the sample data as follows::

        from prolint2 import PL2
        from prolint2.sampledata import COX1

        target_system = PL2(COX1.coordinates, COX1.trajectory)

    """

    def __init__(self):
        super().__init__(
            dataset_name="COX1",
            coordinates_url="https://prolint.github.io/ProLintData/COX1/mid/coordinates.gro",
            trajectory_url="https://prolint.github.io/ProLintData/COX1/mid/trajectory.xtc"
        )


class SMODataSample(_DataSampleBase):
    """
    Class to add sample data in order to test the installation of the package

    It uses the coordinates and trajectory of a simulation of the SMO receptor using the MARTINI force field.

    Example
    -------
    Import the prolint2 library and use the sample data as follows::

        from prolint2 import PL2
        from prolint2.sampledata import SMO

        target_system = PL2(SMO.coordinates, SMO.trajectory)

    """

    def __init__(self):
        super().__init__(
            dataset_name="SMO",
            coordinates_url="https://prolint.github.io/ProLintData/Smoothened/mid/coordinates.gro",
            trajectory_url="https://prolint.github.io/ProLintData/Smoothened/mid/trajectory.xtc"
        )
