#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <degraderthickness_from> <degraderthickness_increment> <degrader_thickness_to>"

if [ $# -ne 3 ]; then
	echo "No arguments supplied, exiting"
	exit
fi

echo degraderthickness_from: $1, degraderthickness_increment: $2 degraderthickness_to: $3

for i in `seq $1 $2 $3`;
do
   sed -i "58s/[2-6]/$i/" GlobalConstants/Constants.h
	root -l -q 'Scripts/findTracking.C'
done
