#!/bin/bash

# This run_* script is a wrapper for L3stats.sh, so it will loop over several
# copes and models. Note that Contrast N for PPI is always PHYS in these models.


# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"


# Change the type of analysis in the " " marks.

for analysis in "act"; do # "act" "ppi_seed-pTPJ"

	# Define the contrast value and the name you would like in the output. 

	analysistype=${analysis}		
		# for copeinfo in "11 offer-unfairness_pmod" "12 social-nonsocial_pmod" "18 phys"; do
		for copeinfo in "1 c_C" "2 c_F" "3 c_S" "4 C_def" "5 C_rec" "6 F_def" "7 F_rec" "8 S_def" "9 S_rec" "10 rec-def" "11 SF-C_face" "12 F-S_rec-def" "13 F-S" "14 F-C" "15 S-C" "16 rec_SocClose" "17 def_SocClose" "18 rec-def_SocClose"; do
		# split copeinfo variable
		set -- $copeinfo
		copenum=$1
		copename=$2

		NCORES=40
		SCRIPTNAME=${maindir}/code/L3stats-trust.sh
		echo $SCRIPTNAME
		while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
			sleep 1s
		done
		bash $SCRIPTNAME $copenum $copename $analysistype &
		echo "$SCRIPTNAME $copenum $copename $analysistype"

	done
done
