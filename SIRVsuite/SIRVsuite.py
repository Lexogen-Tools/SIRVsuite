import argparse as ap
import os
import sys
import logging
from SIRVsuite import __version__

from SIRVsuite.Pipeline.Coverage.SIRVset_coverage import SIRVsuiteCoverage
from SIRVsuite.Pipeline.Concentration.SIRV_concentration import (
    SIRVsuiteConcentration)
from SIRVsuite.Pipeline.Correlation.ERCC_correlation import ERCCcorrelation
from SIRVsuite.Pipeline.helper import *

# set logging
logging.root.setLevel(logging.INFO)
log = logging.getLogger()
hndl = log.handlers[0]
formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
hndl.setFormatter(formatter)

def main():
    modules_to_execute = []

    parser = ap.ArgumentParser()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-i', '--sample-sheet', action='store', required=True, nargs=1,
                               help="specify path to the sample sheet")
    required_args.add_argument('-o', '--output-dir', action='store', required=True, nargs=1,
                               help="specify output directory")
    parser.add_argument('-a', '--all-modules', action='store_true', required=False,
                        help="""this options has the same functionality as
                         --ERCC-correlation --SIRV-concentration --coverage and executes all 3 modules""")
    parser.add_argument('--ERCC-correlation', action='store_true',
                        help="""calculate correlation between 
                        expected and measured concentration for ERCC genes""")
    parser.add_argument('--SIRV-concentration', action='store_true',
                        help="""create heatmap and boxplots
                        for SIRV E0 relative concentration""")
    parser.add_argument('--coverage', action='store_true',
                        help="create coverages for all spike-in genes")
    parser.add_argument('--experiment-name', action='store', default="", required=False, nargs=1,
                        help="specify name of an experiment")
    parser.add_argument('--verbose', action='store_true', required=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # determine which modules should run
    if (args.all_modules or args.SIRV_concentration):
        modules_to_execute.append("SIRV-concentration")
    if (args.all_modules or args.coverage):
        modules_to_execute.append("coverage")
    if (args.all_modules or args.ERCC_correlation):
        modules_to_execute.append("ERCC-correlation")

    if (len(modules_to_execute) == 0):
        raise ValueError("a module has to be specified.. use -h or --help arg for more information")

    if (args.verbose):
        logging.root.setLevel(logging.DEBUG)

    out = os.path.expanduser(args.output_dir[0])
    input_path = os.path.expanduser(args.sample_sheet[0])
    # read sample sheet + module availibility check based on sample sheet
    input_dict = read_sample_sheet(input_path, modules_to_execute=modules_to_execute)

    # run modules
    if "ERCC-correlation" in modules_to_execute or "SIRV-concentration" in modules_to_execute:

        if (args.all_modules or args.SIRV_concentration):
            module_concentration = SIRVsuiteConcentration(sample_sheet=input_dict["SIRV-concentration"], output_dir=out)
            module_concentration.create_sirvsuite_boxplot(module_concentration.data)
            module_concentration.create_sirvsuite_heatmap(module_concentration.data)

        if (args.all_modules or args.ERCC_correlation):
            if (args.SIRV_concentration):
                # If counts have been loaded already for previous module, save computation time
                cnts = module_concentration.cnts
                module_correlation = ERCCcorrelation()
                module_correlation.ERCC_correlation(cnts, output_dir=os.path.join(out, "correlation/"))

                del module_concentration
                del cnts
            else:
                module_correlation = ERCCcorrelation(sample_sheet=input_dict["ERCC-correlation"], output_dir=out)
            
            del module_correlation

    if "coverage" in modules_to_execute:
        module_coverage = SIRVsuiteCoverage(sample_sheet=input_dict["coverage"], output_dir=out, experiment_name=args.experiment_name)
        module_coverage.expected_coverage()
        module_coverage.bam_to_coverage()
        module_coverage.calc_statistics()
        module_coverage.coverage_plot()

        del module_coverage


if __name__ == '__main__':
    sys.exit(main())
