#!/bin/bash
# create date: 2017.02.17, by Nova Smedley
#
# Run the entire pipeline for a given set of cSPADE parameters,
# from SPM pattern generation to logit modeling of patterns.
# There is the option to skip SPM pattern generation, assuming
# the patterns already exist and ready for logit training.
#
# Program executes for each tumor volume type.
# Assumes outDir exists.
#
# Must run in directory of 'run' R scripts
# E.g., '$ ./run_logit.sh --spm no -s 0.3 -g 60 -l 2 -z 2 -d default -p logit_test --dir ~/Data/GBM --saveLogit no -t rate'


while [[ $# -gt 1 ]]
do
key="$1"

case $key in
	-t|--tType) # tumor volume type: rate, c, p, r, or all
	type="$2"
	shift	
	;;	
	--spm) # run spm? yes or no, if yes, need cSPADe constraints
	spm="$2"
	shift
	;;
	-s|--support) # cSPADE support
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
	-d|--dataFolder) # name of data folder (in dir) to read in spm patterns from; 'default' for auto naming
	d="$2"
	shift
	;;	
	-p|--prefix) # name of experiment (where to store logit in outDir)
	p="$2"
	shift
	;;
	--dir) # directory to read spm data and store logit data
	o="$2"
	shift
	;;
	--saveLogit) # save logit data? yes or no, note: results will be saved in a log file regardless
	save="$2"
	shift
	;;
	*)
		# unknown option
	;;
esac
shift # past argument or value
done

# default experiment name
if [[ $d == "default" ]]; then
 d="sup"$s'g'$g'l'$l'z'$z
fi


if [[ $spm == 'yes' ]]; then
	./run_spm.sh --minSupp $s --maxgap $g --maxlength $l --maxsize $z --outDir $o/$p --minSuppList no --tType $type
fi

if [[ $type == 'all' ]]; then
 Rscript --max-ppsize=70000 run_logits_spm.R --tType rate --maxgap $g --maxlength $l --dataFolder $d --prefix $p --dir $o --saveLogit $save
 Rscript --max-ppsize=70000 run_logits_spm.R --tType c --maxgap $g --maxlength $l --dataFolder $d --prefix $p --dir $o --saveLogit $save
 Rscript --max-ppsize=70000 run_logits_spm.R --tType p --maxgap $g --maxlength $l --dataFolder $d --prefix $p --dir $o --saveLogit $save
 Rscript --max-ppsize=70000 run_logits_spm.R --tType r --maxgap $g --maxlength $l --dataFolder $d --prefix $p --dir $o --saveLogit $save
else
 Rscript --max-ppsize=70000 run_logits_spm.R --tType $type --maxgap $g --maxlength $l --dataFolder $d --prefix $p --dir $o --saveLogit $save
fi

