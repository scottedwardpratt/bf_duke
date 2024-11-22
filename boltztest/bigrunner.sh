#! /bin/bash
case $# in
0|1|2)
	echo "Usage: bigrunner.sh run_number first_subrun last_subrun // runs from i=ifirst to <=ifinal" >> crap.txt;
	exit 1 ;;
3)
	cd /home/scott/git/bf_duke/boltztest
	run_number=$1
	rm -f logfiles/run${run_number}*.txt
	first_subrun=$2
	last_subrun=$3
	for((isub=${first_subrun};isub<=${last_subrun};isub++)) do
		runner.sh ${run_number} ${isub}
	done

esac