#!/bin/bash
CORES=50
TRIALS=200
COUNTER=0

BOOTSTRAP=1
BOT_TOP=0
P_OPT=0
MAX_EN=95

while [ ${COUNTER} -lt ${CORES} ];
do

	nohup nice -n 0 /home/efrain/gradient_fit/single_job.sh ${COUNTER} ${TRIALS} ${BOOTSTRAP} ${BOT_TOP} ${P_OPT} ${MAX_EN} &


	let COUNTER+=1
done
