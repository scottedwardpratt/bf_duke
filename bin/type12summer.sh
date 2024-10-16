#!/bin/bash
case $# in
0|1|2)
	echo "Usage: type12summer.sh I J qualifier  // sums for modelruns/runI through modelruns/runJ";
	exit 1 ;;
3)
					
	firsti=$1
	lasti=$2
	qualifier=$3
	#\mkdir -p modelruns/run${ii}/${qualifier}/results_type1_sum
	#\mkdir -p modelruns/run${ii}/${qualifier}/results_type2_sum
	#\mkdir -p modelruns/run${ii}/${qualifier}/results_sum
	for ((ii=${firsti};ii<=${lasti};ii++))
	do
		mkdir -p modelruns/run${ii}/${qualifier}/results_type1_sum
		for jj in `ls modelruns/run${ii}/${qualifier}/results_type1/subruns/`
		do
			for dirname in KK Kp allcharges allcharges_phi0 allcharges_phi45 allcharges_phi90 piK pip pipi pp
			do
				echo ${jj}
				cp -f modelruns/run${ii}/${qualifier}/results_type1/subruns/${jj}/${dirname}
			done
		done
										
	done
esac