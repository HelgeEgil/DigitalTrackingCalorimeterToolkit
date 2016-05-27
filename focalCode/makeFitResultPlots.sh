#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <energy_from> <energy_increment> <energy_to>"

if [ $# -ne 3 ]; then
	echo "No arguments supplied, exiting"
	exit
fi

echo Energy_from: $1, energy_increment: $2 energy_to: $3

if [ $1 -gt $3 ] ; then
	echo "energy_from ($2) is higher than energy_to ($4)"
	echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to>"
	exit
fi

for i in `seq $1 $2 $3`;
do
	echo "Making plots at $i MeV."
	# makeTracksFile.C(int Runs, int dataType, int energy)
	root -l -q 'Scripts/makeBraggPeakPDF.C('$i')'
done
