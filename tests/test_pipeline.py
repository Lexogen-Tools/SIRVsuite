import pytest
from SIRVsuite.Pipeline.helper import read_sample_sheet


'''
@pytest.fixture
def empty_sample_sheet():
    return "./test_data/sample_sheet_fail/sample_sheet_empty.tsv"

@pytest.fixture
def wrong_sample_sheet():
    return "./test_data/sample_sheet_fail/sample_sheet_wrong.tsv"

def test_sample_sheet(empty_sample_sheet):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        read_sample_sheet(empty_sample_sheet)
    assert pytest_wrapped_e.type == SystemExit
'''

@pytest.mark.parametrize("sample_sheet_pass", [
    "./test_data/sample_sheet_pass/sample_sheet_test_SIRVset3.tsv",
    "./test_data/sample_sheet_pass/sample_sheet_test_SIRVset4.tsv"
])

def test_sample_pass(sample_sheet_pass):
	read_sample_sheet(sample_sheet_pass)
