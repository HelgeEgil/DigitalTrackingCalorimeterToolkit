#!/bin/bash

NCORES=4
IDX=1

for j in `seq -w 35 1 35`; do
   for i in `seq -w 210 1 210`; do
      hdt=`echo "scale=3; -$i/2-2" | bc`
      beampos=`echo "scale=3; -$i-5" | bc`
      Gate -a "'[absorberthickness,$j] [energy,250] [degraderthickness,$i] [halfdegraderthickness,$hdt] [beampos,$beampos]'" Main.mac > terminal_output.txt &
      PIDLIST="$PIDLIST $!"
      echo "Running: Absorber thickness = $j mm; degrader thickness = $i mm"

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
