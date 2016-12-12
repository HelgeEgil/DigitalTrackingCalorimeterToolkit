#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

echo "Usage: ./run.sh <energy_from> <energy_step> <energy_to> <npart>"
echo Energy From: $1, Energy Step: $2, Energy To: $3, npart: $4

if [ $# -ne 4 ]; then
	echo Invalid number of arguments: $#
	exit
fi
for i in `seq $1 $2 $3`;
do
   nice -n 10 Gate -a "'[npart,$4] [energy,$i]'" Main.mac > terminal_output.txt &
   echo "Running: $i"
done
