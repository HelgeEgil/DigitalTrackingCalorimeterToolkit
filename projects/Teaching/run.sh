#!/bin/bash

echo "Usage: ./run.sh <energy_from> <energy_increment> <energy_to>"
echo Energy_from: $1, energy_increment: $2, energy_to: $3

if [ $# -ne 3 ]; then
	echo Invalid number of arguments: $#
	exit
fi

for i in `seq $1 $2 $3`;
do
	nice -n 10 Gate -a "'[energy,$i]" waterphantom.mac > terminal_output.txt &
	echo "Running: $i"
done
