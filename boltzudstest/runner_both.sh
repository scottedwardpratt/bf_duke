#!/bin/bash
case $# in
	0|1)
		echo "Usage: runner_both.sh ifirst ifinal  // runs from i=ifirst to <=ifinal";
  	exit 1 ;;
	2)
		firsti=$1
		lasti=$2
		for ((ii=${firsti};ii<=${lasti};ii++))
		do
			echo "____________ hydro2uds/b3d for run number " ${ii} ______________;
			../bin/hydro2uds $ii;
			../bin/b3d default_${ii} ${ii} ${ii};
		done
esac
