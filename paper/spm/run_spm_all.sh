#!/bin/bash
# create date: 2017.02.17, by Nova Smedley
#
# Calls 'run_spm.sh'.
# Call it multiple times with different cSPADE 
# parameters (predefined here) for a given maxlength argument.
#
# Must run in directory of 'run' R scripts
# E.g., '$ ./run_spm_all.sh -l 2 -t all -o ~/Data/GBM'
#
# Assumes outDir exists.

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
	-l|--maxlength) # cSPADE
	l="$2"
	shift
	;;
	-t|--tType) # tumor volume type: rate, c, p, r, or all
	type="$2"
	shift
	;;	
	-o|--outDir) # directory to store exerpiment data
	o="$2"
	shift
	;;
	*)
		# unknown option
	;;
esac
shift # past argument or value
done


minSupp=0.2

for z in 2 3 4  # there are only 7 categories of clinical events
do
	./run_spm.sh -s $minSupp --tType $type -g 60 -l $l -z $z -o $o
	./run_spm.sh -s $minSupp --tType $type -g 45 -l $l -z $z -o $o
	./run_spm.sh -s $minSupp --tType $type -g 30 -l $l -z $z -o $o
	./run_spm.sh -s $minSupp --tType $type -g 15 -l $l -z $z -o $o
done


