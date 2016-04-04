#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to> <npart> <sigma> <thickness>"
echo Energy_from: $1, energy_increment: $2, energy_to: $3, npart: $4, accuracy um: $5, sigma MeV: $6, thickness: $7

if [ $# -ne 7 ]; then
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
	time Gate -a "'[energy,$i] [npart,$4] [step_active,$5] [sigma,$6] [thick,$7]'" focal.mac > terminal_output.txt &
done
