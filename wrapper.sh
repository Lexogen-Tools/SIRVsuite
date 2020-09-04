#/bin/bash

usage() { echo "Usage:" $0 [-i <sample_sheet>] [-o output_dir] [-] 1>&2 exit 1; }
 
while getopts ":i:o:" o; do
	case "${o}" in 
		i)
			input=${OPTARG}
			echo "input: ${input}"
			;;
		o) 
			output=${OPTARG}
			echo "output: ${output}"
			;;
		*) usage
			;;
	esac
done
