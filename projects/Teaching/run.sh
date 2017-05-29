#!/bin/bash

echo "Usage: ./run.sh <degraderthickness_from> <degraderthickness_increment> <degraderthickness_to>"
echo Absorber thickness: $1, Phantom thickness from $2 step $3 to $4

NCORES=10
IDX=1

if [ $# -ne 3 ]; then
	echo Invalid number of arguments: $#
	exit
fi

for i in `seq $1 $2 $3`; do
	nice -n 10 Gate -a "'[energy,$i]'" waterphantom.mac > terminal_output.txt &
   echo "Running: $i"
   if (( $IDX % $NCORES == 0 )); then 
      echo "Waiting for run $i (PID $!)"
      IDX=1
      time wait $!
   else
      IDX=$(( IDX+1 ))
   fi
done
