#!/bin/bash

NCORES=4
IDX=1

for x in `seq -w 15 5 15`; do
   for y in `seq -w 15 5 15`; do
      Gate -a "'[absorberthickness,3] [energy,1800] [spotx,$x] [spoty,$y]'" Main_Head.mac > terminal_output.txt &
#      Gate -a "'[absorberthickness,3] [energy,600] [spotx,$x] [spoty,$y]'" Main_Head.mac > terminal_output.txt &
      PIDLIST="$PIDLIST $!"
      echo "Running: Spot position ($x,$y)"

      if (( $IDX % $NCORES == 0 )); then 
         echo "Waiting for (PIDS $PIDLIST)"
         time wait $PIDLIST
         unset PIDLIST
         IDX=1
      else
         IDX=$(( IDX+1 ))
      fi
   done

   echo "Waiting for (PIDS $PIDLIST)"
   time wait $PIDLIST
   unset PIDLIST
   IDX=1
done

echo "Waiting for (PIDS $PIDLIST)"
time wait $PIDLIST
unset PIDLIST
IDX=1
