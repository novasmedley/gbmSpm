#!/bin/bash
# create date: 2017.02.17, by Nova Smedley
#
#
# Calls 'runSPM.R'.
# Call it multiple times, once for each tumor volume types 
# predefined in 'run_spm.R' for a given set of cSPADE parameters.
#
# Note that 'run_spm.R' will take minSupp and by default will
# also generate a sequence of min supports with the remaining cSPADE parameters.
# Therefore cSPADE is executed multiple times for each min support.
# This can be turned off by the argument '--minSuppList' with 'no.'
#
# Must run in directory of 'run' R scripts
# E.g., '$ ./run_spm.sh -o ~/Data/GBM -s 0.4 -g 60 -l 2 -z 2 -t rate'
#
# Assumes outDir exists.

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
	-t|--tType) # tumor volume type: rate, c, p, r, or all
	type="$2"
	shift
	;;	
	-s|--minSupp) # cSPADE
	s="$2"
	shift
	;;
	-g|--maxgap) # cSPADE
	g="$2"
	shift
	;;
	-l|--maxlength) # cSPADE
	l="$2"
	shift
	;;
	-z|--maxsize) # cSPADE
	z="$2"
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

if [[ $type == "all" ]]; then
 Rscript --max-ppsize=70000 run_spm.R --tType rate --minSupp $s --maxgap $g --maxlength $l --maxsize $z --outDir $o --minSuppList "yes"
 Rscript --max-ppsize=70000 run_spm.R --tType c --minSupp $s --maxgap $g --maxlength $l --maxsize $z --outDir $o --minSuppList "yes"
 Rscript --max-ppsize=70000 run_spm.R --tType p --minSupp $s --maxgap $g --maxlength $l --maxsize $z --outDir $o	--minSuppList "yes"
 Rscript --max-ppsize=70000 run_spm.R --tType r --minSupp $s --maxgap $g --maxlength $l --maxsize $z --outDir $o	--minSuppList "yes"
else
 Rscript --max-ppsize=70000 run_spm.R --tType $type --minSupp $s --maxgap $g --maxlength $l --maxsize $z --outDir $o --minSuppList "yes"
fi

