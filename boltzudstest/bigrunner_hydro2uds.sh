#!/bin/bash
nproc=24
nruns=30
for ((i=0;i<${nproc};i+=1))
do
	firsti=`expr ${i} \* ${nruns}`;
	lasti=`expr ${firsti} + ${nruns} - 1`;
	rm -f logfiles/hydro2uds_${firsti}_${lasti}.txt;
	echo starting runs with firsti=${firsti}, lasti=${lasti};
	`./runner_hydro2uds.sh ${firsti} ${lasti} > logfiles/hydro2uds_${firsti}_${lasti}.txt &` ;
done
