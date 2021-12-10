from _pytest.fixtures import fixture
import pytest
from SIRVsuite.Pipeline.helper import read_sample_sheet
from SIRVsuite.Pipeline.countReader import countReader
from SIRVsuite.Pipeline.Coverage.SIRVset_coverage import SIRVsuiteCoverage as cov_module
import os
import shutil

# Check intented PASS when providing wrong transcript count table
@pytest.mark.parametrize("count_tables_pass", 
[
    "./tests/test_data/mix2/SIRVset3/test_transcript_count_1",
    "./tests/test_data/mix2/SIRVset4/test_transcript_count_4",
    "./tests/test_data/mix2/test_tcount_shuffle_col"
])

def test_count_read_pass(count_tables_pass):
    reader = countReader()
    cnts = reader.read_counting_file([count_tables_pass], counting_method="mix2")
    # Check if no spike-in transript is missing (but take into consideration missing SIRV108 for certain Lot. numbers)
    assert (len(cnts["ERCC"]), len(cnts["SIRV"])) in [(92, 68), (92, 69)]

def test_count_htseq():
    reader = countReader()
    cnts = reader.read_counting_file(["./tests/test_data/mix2/test_tcount_htseq"], counting_method="htseq")
    assert (len(cnts["ERCC"]), len(cnts["SIRV"])) in [(92, 68), (92, 69)]

# Check intented FAIL when providing wrong transcript count table
@pytest.mark.parametrize("count_tables_fail", 
[
    "./tests/test_data/mix2/test_tcount_fail",
])

def test_count_read_fail(count_tables_fail):
    with pytest.raises(ValueError) as pytest_wrapped_e:
        reader = countReader()
        _ = reader.read_counting_file([count_tables_fail])
    assert pytest_wrapped_e.type == ValueError

# Test the whole pipeline

@pytest.mark.parametrize("module_params", 
[
    "--ERCC-correlation",
    "--SIRV-concentration",
    "--ERCC-correlation --SIRV-concentration",
    "--coverage"
])


def test_entrypoint(module_params):
    exit_status = os.system("PYTHONPATH=../ python -m SIRVsuite.SIRVsuite -i ./tests/test_data/sample_sheet_pass/sample_sheet_test_SIRVset3.tsv -o ./tests/test_data/temp/ %s"%(module_params))
    if os.path.exists("./tests/test_data/temp"):
        shutil.rmtree("./tests/test_data/temp")
    assert exit_status == 0
