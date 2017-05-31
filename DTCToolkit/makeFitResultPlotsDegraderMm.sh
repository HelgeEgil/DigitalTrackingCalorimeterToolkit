#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <degraderthickness_from> <degraderthickness_increment> <degrader_thickness_to> <absorberthickness_from> <absorberthickness_increment> <absorberthickness_to>"

if [ $# -ne 6 ]; then
	echo "No arguments supplied, exiting"
	exit
fi

echo degraderthickness_from: $1, degraderthickness_increment: $2 degraderthickness_to: $3
echo absorberthickness_from: $4, absorberthickness_increment: $5 absorberthickness_to: $6

for j in `seq $4 $5 $6`;
do
   sed -i "58s/[2-6]/$j/" GlobalConstants/Constants.h
   for i in `seq $1 $2 $3`;
   do
   	echo "Making plots at $i mm degrader with $j mm absorber."
   	# makeTracksFile.C(int Runs, int dataType, int energy)
	   #root -l -q 'Scripts/makeBraggPeakPDFDegrader.C('$i')'
   done
done
