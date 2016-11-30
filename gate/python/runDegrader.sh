#!/bin/bash

echo "Usage: ./run.sh <absorberthickness> <degraderthickness_from> <degraderthickness_increment> <degraderthickness_to>"
echo Absorber thickness: $1, degraderthickness from $2 step $3 to $4

if [ $# -ne 4 ]; then
	echo Invalid number of arguments: $#
	exit
fi

for i in `seq $2 $3 $4`;
do
	nice -n 10 Gate -a "'[absorberthickness,$1] [energy,230] [degraderthickness,$i]" Main.mac > terminal_output.txt &
	echo "Running: $i"
done
