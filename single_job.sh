#!/bin/bash
#nohup nice -n 0 ./home/efrain/gradient_fit/build/single_job.sh ${COUNTER} ${TRIALS} ${BOOTSTRAP} ${BOT_TOP} ${P_OPT} ${MAX_EN} &


for (( RUN=0; RUN<$2; RUN++ ))
do
	THIS_TRIAL=$(($1*($2)+RUN))

	/home/efrain/gradient_fit/build/gsl_trial $3 $4 $5 $6 > /xdata/depolarization_fits/$3_$4_$5_$6/${THIS_TRIAL}.txt 

done
