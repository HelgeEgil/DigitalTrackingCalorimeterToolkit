#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <material> <energy_from> <energy_increment> <energy_to> <npart> <sigma_energy>"
echo Material: $1, energy_from: $2, energy_increment: $3, energy_to: $4, npart: $5, sigma_energy: $6

if [ $# -ne 6 ]; then
	echo Invalid number of arguments: $#
	exit
fi

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

for i in `seq $2 $3 $4`;
do
	nice -n 15 Gate -a "'[material,$1] [energy,$i] [npart,$5] [step_active,10] [step_passive,80] [sigma,$6]'" focal.mac > terminal_output.txt &
done
