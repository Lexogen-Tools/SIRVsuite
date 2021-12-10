from _pytest.fixtures import fixture
import pytest

from SIRVsuite.Pipeline.Coverage.SIRVset_coverage import SIRVsuiteCoverage as cov_module
import os
import shutil

def test_coverage_module():
    
    covm = cov_module(
        sample_sheet={'test_sample': {
            "alignment_path": "./tests/test_data/alignment/SIRVset3/test_alignment_1.bam",
            "read_orientation": "fwd"
        }},
        output_dir="./tests/test_data/temp/test_cov1", experiment_name="")

    # Check if stat module returns an error when called without expected coverage
    with pytest.raises(Exception):
        covm.calc_statistics()

    # Check for errors for regular expected coverage and with transient lengths
    covm.expected_coverage(transition_lengths=True)
    covm.expected_coverage()

    # Check if stat module returns an error when called without bam coverage
    with pytest.raises(Exception):
        covm.calc_statistics()

    # Check for error while calculating bam coverage and exporting it to COV or BW
    covm.bam_to_coverage(output_type="cov")
    covm.bam_to_coverage(output_type="BigWig")

    # Check for errors in statistics module
    covm.calc_statistics()

    # Check for errors during coverage plot generation
    covm.coverage_plot()

    covm = cov_module(
        sample_sheet={'test_sample': {
            "alignment_path": "./tests/test_data/alignment/SIRVset3/test_alignment_1.bam",
            "read_orientation": "fwd"
        }},
        output_dir="./tests/test_data/temp/test_cov2", experiment_name="test experiment")
    
    covm.expected_coverage()
    covm.bam_to_coverage()
    covm.calc_statistics()
    covm.coverage_plot()

    # Check for erros when converting annotation to UTR annotation
    covm.annot_2_UTR()

    with pytest.raises(Exception):
        covm.load_annotation("./tests/non_existing_file.gtf")

    if os.path.exists("./tests/test_data/temp"):
       shutil.rmtree("./tests/test_data/temp")