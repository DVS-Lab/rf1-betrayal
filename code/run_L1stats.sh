#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"
nruns=2
task=ugr # edit if necessary

for ppi in "pTPJ"; do #0
#for ppi in "ecn"; do # putting 0 first will indicate "activation" Put in "NAcc-bin" and "ecn" for PPI 
#for ppi in "pTPJ"; do
	 for sub in 10317; do
	 #for sub in `cat ${scriptdir}/sublist_allUGRresponses.txt`; do	
	#for sub in `cat ${scriptdir}/sublist_mkExist.txt`; do	
	  for run in `seq $nruns`; do

	  	# Manages the number of jobs and cores
	  	SCRIPTNAME=${basedir}/code/L1stats.sh
	  	NCORES=30
	  	while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
	    		sleep 5s
	  	done
	  	bash $SCRIPTNAME $sub $run $ppi &
			sleep 1s
	  done
	done
done
