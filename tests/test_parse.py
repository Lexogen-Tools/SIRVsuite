import pytest
import glob
from SIRVsuite.Pipeline.helper import read_sample_sheet

"""
@pytest.mark.parametrize("sample_sheet_fail", [
    "./tests/test_data/sample_sheet_fail/sample_sheet_empty.tsv",
    "./tests/test_data/sample_sheet_fail/sample_sheet_wrong.tsv",
])
"""

@pytest.mark.parametrize("sample_sheet_fail", glob.glob("./tests/test_data/sample_sheet_fail/*.tsv"))

def test_sample_fail(sample_sheet_fail):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        read_sample_sheet(sample_sheet_fail)
    assert pytest_wrapped_e.type == SystemExit

@pytest.mark.parametrize("sample_sheet_pass", glob.glob("./tests/test_data/sample_sheet_pass/*.tsv"))

def test_sample_pass(sample_sheet_pass):
	read_sample_sheet(sample_sheet_pass)
