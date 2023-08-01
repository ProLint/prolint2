r"""UFCC sample data
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os
import requests
import tqdm


def file_download(url: str, filename: str, desc: str):
    with open(filename, "wb") as f:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))

            tqdm_params = {
                "desc": url,
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


class GIRKDataSample:
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
        self.path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "data/GIRK/"
        )
        if (
            os.path.isfile(os.path.join(self.path, "coordinates.gro")) is False
            or os.path.isfile(os.path.join(self.path, "trajectory.xtc")) is False
        ):
            file_download(
                "https://prolint.github.io/ProLintData/GIRK/mid/coordinates.gro",
                os.path.join(self.path, "coordinates.gro"),
                "Downloading GIRK coordinates...",
            )
            file_download(
                "https://prolint.github.io/ProLintData/GIRK/mid/trajectory.xtc",
                os.path.join(self.path, "trajectory.xtc"),
                "Downloading GIRK trajectory...",
            )
            self.coordinates = os.path.join(self.path, "coordinates.gro")
            self.trajectory = os.path.join(self.path, "trajectory.xtc")
        else:
            self.coordinates = os.path.join(self.path, "coordinates.gro")
            self.trajectory = os.path.join(self.path, "trajectory.xtc")
        self.contacts = os.path.join(self.path, "contacts.csv")


class COX1DataSample:
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
        self.path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "data/COX1/"
        )
        if (
            os.path.isfile(os.path.join(self.path, "coordinates.gro")) is False
            or os.path.isfile(os.path.join(self.path, "trajectory.xtc")) is False
        ):
            file_download(
                "https://prolint.github.io/ProLintData/COX1/mid/coordinates.gro",
                os.path.join(self.path, "coordinates.gro"),
                "Downloading COX1 coordinates...",
            )
            file_download(
                "https://prolint.github.io/ProLintData/COX1/mid/trajectory.xtc",
                os.path.join(self.path, "trajectory.xtc"),
                "Downloading COX1 trajectory...",
            )
            self.coordinates = os.path.join(self.path, "coordinates.gro")
            self.trajectory = os.path.join(self.path, "trajectory.xtc")
        else:
            self.coordinates = os.path.join(self.path, "coordinates.gro")
            self.trajectory = os.path.join(self.path, "trajectory.xtc")
        self.contacts = os.path.join(self.path, "contacts.csv")


class SMODataSample:
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
        self.path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "data/SMO/"
        )
        if (
            os.path.isfile(os.path.join(self.path, "coordinates.gro")) is False
            or os.path.isfile(os.path.join(self.path, "trajectory.xtc")) is False
        ):
            file_download(
                "https://prolint.github.io/ProLintData/Smoothened/mid/coordinates.gro",
                os.path.join(self.path, "coordinates.gro"),
                "Downloading SMO coordinates...",
            )
            file_download(
                "https://prolint.github.io/ProLintData/Smoothened/mid/trajectory.xtc",
                os.path.join(self.path, "trajectory.xtc"),
                "Downloading SMO trajectory...",
            )
            self.coordinates = os.path.join(self.path, "coordinates.gro")
            self.trajectory = os.path.join(self.path, "trajectory.xtc")
        else:
            self.coordinates = os.path.join(self.path, "coordinates.gro")
            self.trajectory = os.path.join(self.path, "trajectory.xtc")
        self.contacts = os.path.join(self.path, "contacts.csv")
