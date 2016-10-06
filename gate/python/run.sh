#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <absorberthickness> <energy_from> <energy_increment> <energy_to>"
echo Material: $1, energy_from: $2, energy_increment: $3, energy_to: $4, npart: $5, sigma_energy: $6

if [ $# -ne 4 ]; then
	echo Invalid number of arguments: $#
	exit
fi

for i in `seq $2 $3 $4`;
do
	time nice -n 10 Gate -a "'[absorberthickness,$1] [energy,$i]" Main.mac > terminal_output.txt &
	echo "Running: $i"
done
