#!/bin/bash

echo "Usage: ./run.sh <energy_from> <energy_increment> <energy_to>"
echo Energy from $1 step $2 to $3

NCORES=10
IDX=1

if [ $# -ne 3 ]; then
	echo Invalid number of arguments: $#
	exit
fi

for i in `seq $1 $2 $3`; do
	time nice -n 10 Gate -a "'[energy,$i]'" waterphantom.mac > terminal_output.txt &
   echo "Running: $i"
   if (( $IDX % $NCORES == 0 )); then 
      echo "Waiting for run $i (PID $!)"
      IDX=1
      time wait $!
   else
      IDX=$(( IDX+1 ))
   fi
done
