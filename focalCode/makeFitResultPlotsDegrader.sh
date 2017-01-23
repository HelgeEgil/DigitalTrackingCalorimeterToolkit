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
	echo "Making plots at $i mm."
	# makeTracksFile.C(int Runs, int dataType, int energy)
	root -l -q 'Scripts/makeBraggPeakPDFDegrader.C('$i')'
done
