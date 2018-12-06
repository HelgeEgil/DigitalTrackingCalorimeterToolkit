#!/bin/bash

NCORES=4
IDX=1

for j in `seq -w 3 1 3`; do
   for i in `seq -w 100 1 103`; do
      hdt=`echo "scale=3; -$i/2-2" | bc`
      beampos=`echo "scale=3; -$i-5" | bc`
      Gate -a "'[absorberthickness,$j] [energy,230] [degraderthickness,$i] [halfdegraderthickness,$hdt] [beampos,$beampos]'" Main_full.mac > terminal_output.txt &
      PIDLIST="$PIDLIST $!"
      echo "Running (full): Absorber thickness = $j mm; degrader thickness = $i mm"
      if (( $IDX % $NCORES == 0)); then
         echo "Waiting for (PIDS $PIDLIST)"
         time wait $PIDLIST
         unset PIDLIST
         IDX=1
      else
         IDX=$((IDX+1))
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
