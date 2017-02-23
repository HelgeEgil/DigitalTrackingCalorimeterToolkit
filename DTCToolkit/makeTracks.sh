#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <energy_from> <energy_increment> <energy_to> <runs>"

if [ $# -ne 4 ]; then
	echo "No arguments supplied, exiting"
	exit
fi

echo Material: energy_from: $1, energy_increment: $2, energy_to: $3, npart: $4

if [ $1 -gt $3 ] ; then
	echo "energy_from ($1) is higher than energy_to ($3)"
	echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to>"
	exit
fi

for i in `seq $1 $2 $3`;
do
	root -l -q 'Scripts/makeTracksFile.C('$4', '$i')'
done
