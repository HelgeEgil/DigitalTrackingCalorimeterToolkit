#!/bin/bash

NCORES=6
IDX=1

for x in `seq -w -24 6 72`; do
   for y in `seq -w -72 6 72`; do
      Gate -a "'[absorberthickness,3] [energy,2280] [spotx,$x] [spoty,$y]'" Main_Head_Carbon.mac > terminal_output.txt &
      PIDLIST="$PIDLIST $!"
      echo "Running Carbon: Spot position ($x,$y)"

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

for x in `seq -w -24 6 72`; do
   for y in `seq -w -30 6 72`; do
      Gate -a "'[absorberthickness,3] [energy,760] [spotx,$x] [spoty,$y]'" Main_Head_Helium.mac > terminal_output.txt &
      PIDLIST="$PIDLIST $!"
      echo "Running Helium: Spot position ($x,$y)"

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
