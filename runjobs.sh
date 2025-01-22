#!/bin/bash
CORES=50
TRIALS=200
COUNTER=0

BOOTSTRAP=1
BOT_TOP=0
P_OPT=0
MAX_EN=165
DATASET=1

while [ ${COUNTER} -lt ${CORES} ];
do

	nohup nice -n 0 /home/efrain/gradient_fit/single_job.sh ${COUNTER} ${TRIALS} ${BOOTSTRAP} ${BOT_TOP} ${P_OPT} ${MAX_EN} ${DATASET} &


	let COUNTER+=1
done



#	Wrong number of arguments used.
#	        Please instead use: ./code [bootstrap opt] [bot_top opt] [diffuse opt] [emax]
#	**********************************
#	                [bootstrap opt == 0]: Don't bootstrap
#	                [bootstrap opt == 1]: Bootstrap
#	**********************************
#	                [bot_top opt == 0]: Fit bottom chamber
#	                [bot_top opt == 1]: Fit top chamber
#	**********************************
#	                [diffuse opt == 0]: p_ring=p_electrode
#	                [diffuse opt == 1]: p_ring=0
#	                [diffuse opt == 2]: p_ring=1
#	                [diffuse opt == 3]: p_ring, p_electrode
#	**********************************
#	                [emax]: 35 (rexolite) or 54 (aluminium)
#	**********************************
#			[rexolite quartz opt == 0]: Fit rexolite data
#			[rexolite quartz opt == 1]: Fit quartz data
#	**********************************
