#!/bin/bash
case $# in
	0|1)
		echo "Usage: runner_hydro2uds.sh ifirst ifinal  // runs from i=ifirst to <=ifinal";
  	exit 1 ;;
	2)
		firsti=$1
		lasti=$2
		for ((ii=${firsti};ii<=${lasti};ii++))
		do
			echo "____________ hydro2uds for run number " ${ii} ______________;
			./hydro2uds $ii;
		done
esac