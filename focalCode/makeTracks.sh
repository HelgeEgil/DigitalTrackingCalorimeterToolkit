#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to> <npart>"

if [ $# -ne 5 ]; then
	echo "No arguments supplied, exiting"
	exit
fi

echo Material: $1, energy_from: $2, energy_increment: $3, energy_to: $4, npart: $5

if [ "$1" != "Aluminium" ] && [ "$1" != "Tungsten" ] && [ "$1" != "myPMMA" ] ; then
	echo Please input a valid material. Your choice: $1
	echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to>"
	exit
fi

if [ $2 -gt $4 ] ; then
	echo "energy_from ($2) is higher than energy_to ($4)"
	echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to>"
	exit
fi

for i in `seq --format="%.2f" $2 $3 $4`;
do
	echo "Moving .root file"
	mv rawdata/test_"$1"_"$i"MeV.root rawdata/test.root
	echo "Making tracks at $i MeV."
	# makeTracksFile.C(int Runs, int dataType, int energy)
	root -l -q 'makeTracksFile.C('$5', '$i')'
	mv rawdata/test.root rawdata/test_"$1"_"$i"MeV.root
done
