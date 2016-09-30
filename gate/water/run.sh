#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to> <npart> <sigma>"
echo Energy_from: $1, energy_increment: $2, energy_to: $3, npart: $4, accuracy um: $5

if [ $# -ne 6 ]; then
	echo Invalid number of arguments: $#
	exit
fi

if [ $1 -gt $3 ] ; then
	echo "energy_from ($1) is higher than energy_to ($3)"
	echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to>"
	exit
fi

for i in `seq $1 $2 $3`;
do
	Gate -a "'[energy,$i] [npart,$4] [step_active,$5] [sigma,$6]'" focal_script.mac > terminal_output.txt &
done
