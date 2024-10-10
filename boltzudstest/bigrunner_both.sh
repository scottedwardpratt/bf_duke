#! /bin/bash
nproc=24
nruns=30
for ((i=0;i<${nproc};i+=1))
do
	firsti=`expr ${i} \* ${nruns}`;
	lasti=`expr ${firsti} + ${nruns} - 1`;
	rm -f logfiles/both_${firsti}_${lasti}.txt;
	echo starting runs with firsti=${firsti}, lasti=${lasti};
	`./runner_both.sh ${firsti} ${lasti} > logfiles/both_${firsti}_${lasti}.txt &` ;
done
