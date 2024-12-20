#!/bin/bash
case $# in
	0|1)
		echo "Usage: runner.sh run_number subrun_number // runs from i=ifirst to <=ifinal" >> crap.txt;
  	exit 1 ;;
	2)
		run_number=$1
		subrun_number=$2
		rm -f logfiles/run${run_number}_subrun${subrun_number}.txt
		../bin/spectrav2 ${run_number} ${subrun_number} > logfiles/runscript_run${run_number}_subrun${subrun_number}.txt

esac