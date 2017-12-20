#!/bin/bash

echo "Usage: ./run.sh <degraderthickness_from> <degraderthickness_increment> <degraderthickness_to>"
echo Phantom thickness from $1 step $2 to $3

NCORES=8
IDX=1

if [ $# -ne 3 ]; then
	echo Invalid number of arguments: $#
	exit
fi

for i in `seq $1 $2 $3`;
do
   hdt=`echo "scale=3; -$i/2-2" | bc`
   beampos=`echo "scale=3; -$i-5" | bc`
	nice -n 10 Gate -a "'[energy,230] [degraderthickness,$i] [halfdegraderthickness,$hdt] [beampos,$beampos]'" waterphantom.mac > terminal_output.txt &
	echo "Running: $i"
   if (( $IDX % $NCORES == 0 )); then
      echo "Waiting for run $i (PID $!)"
      IDX=1
      time wait $!
   else
      IDX=$((IDX+1))
   fi
done
