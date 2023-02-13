import os
from prolint2.sampledata import GIRKDataSample

def test_GIRKDataSample_path():
    girk = GIRKDataSample()
    assert os.path.isdir(girk.path)

def test_GIRKDataSample_coordinates():
    girk = GIRKDataSample()
    assert os.path.isfile(girk.coordinates)

def test_GIRKDataSample_trajectory():
    girk = GIRKDataSample()
    assert os.path.isfile(girk.trajectory)
