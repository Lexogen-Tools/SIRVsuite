import argparse as ap
import os
import sys
import logging

from Pipeline.Coverage.SIRVsuiteCoverage import SIRVsuiteCoverage
from Pipeline.Concentration.SIRVsuite_concentration import (
    SIRVsuiteConcentration)
from Pipeline.Correlation.ERCC_correlation import ERCCcorrelation
from Pipeline.helper import *

# set logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger()


def main():
    modules_to_execute = []

    parser = ap.ArgumentParser()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-i', '--sample-sheet', action='store', required=True, nargs=1,
                               help="Specify path to the sample sheet.")
    required_args.add_argument('-o', '--output-dir', action='store', required=True, nargs=1,
                               help="Specify output directory.")
    parser.add_argument('-a', '--all-modules', action='store_true', required=False,
                        help="""This options has the same functionality as
                         --ERCC-correlation --SIRV-concentration --coverage and executes all 3 modules.""")
    parser.add_argument('--ERCC-correlation', action='store_true',
                        help="""Specify to calculate correlation between 
                        expected and measured concentration for ERCC genes.""")
    parser.add_argument('--SIRV-concentration', action='store_true',
                        help="""Specify to create heatmap and boxplots
                        for SIRV E0 relative concentration.""")
    parser.add_argument('--coverage', action='store_true',
                        help="Specify to create coverages for all spike-in genes.")
    parser.add_argument('--experiment-name', action='store', default="", required=False, nargs=1,
                        help="Specify name of an experiment.")
    args = parser.parse_args()

    # determine which modules should run
    if (args.all_modules or args.ERCC_correlation or args.SIRV_concentration):
        modules_to_execute.append("concentration")

    if (args.coverage or args.all_modules):
        modules_to_execute.append("coverage")

    if (len(modules_to_execute) == 0):
        log.warning("You have not specified any module..")
        parser.print_help()

    out = os.path.expanduser(args.output_dir[0])
    input_path = os.path.expanduser(args.sample_sheet[0])

    # read sample sheet + module availibility check based on sample sheet
    input_dict = read_sample_sheet(input_path, modules_to_execute=modules_to_execute)

    # run modules
    if "concentration" in modules_to_execute:

        a = []
        b = []
        cnts = []

        if (args.all_module or args.SIRV_concentration):
            b = SIRVsuiteConcentration(sample_sheet=input_dict["concentration"], output_dir=out)
            b.create_sirvsuite_boxplot(b.data)
            b.create_sirvsuite_heatmap(b.data)

        if (args.all_module or args.ERCC_correlation):
            if (args.SIRV_concentration):
                cnts = b.cnts
                a = ERCCcorrelation()
                a.ERCC_correlation(cnts, output_dir=os.path.join(out, "correlation/"))
            else:
                a = ERCCcorrelation(sample_sheet=input_dict["concentration"], output_dir=out)

        # delete after processsing
        del cnts
        del a
        del b

    if "coverage" in modules_to_execute:
        c = SIRVsuiteCoverage(sample_sheet=input_dict["coverage"], output_dir=out, experiment_name="")
        c.expected_coverage()
        c.bam_to_coverage()
        c.calc_statistics()
        c.coverage_plot()

        del c


if __name__ == '__main__':
    sys.exit(main())
