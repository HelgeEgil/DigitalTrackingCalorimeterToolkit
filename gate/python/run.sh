#!/bin/bash

echo "Usage: ./run.sh <absorberthickness> <energy_from> <energy_increment> <energy_to>"
echo Absorber thickness: $1, energy_from: $2, energy_increment: $3, energy_to: $4

if [ $# -ne 4 ]; then
	echo Invalid number of arguments: $#
	exit
fi

for i in `seq $2 $3 $4`;
do
	nice -n 10 Gate -a "'[absorberthickness,$1] [energy,$i]" Main.mac > terminal_output.txt &
	echo "Running: $i"
done
