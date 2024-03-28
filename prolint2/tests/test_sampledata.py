import os
import pytest
from prolint2.sampledata import GIRKDataSample, COX1DataSample, SMODataSample

GIRK = GIRKDataSample()
COX1 = COX1DataSample()
SMO = SMODataSample()


@pytest.mark.parametrize("DataSample", [GIRK, COX1, SMO])
def test_sampledata(DataSample):
    assert os.path.isdir(DataSample.path)
    assert os.path.isfile(DataSample.coordinates)
    assert os.path.isfile(DataSample.trajectory)
    assert os.path.isfile(DataSample.contacts)
