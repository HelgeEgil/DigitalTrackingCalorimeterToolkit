#!/bin/bash

NCORES=4
IDX=1

for i in `seq -w 0.2 0.4 15`; do
   ztracker=`echo "scale=3; 160/2+15" | bc`
   zbeam=`echo "scale=3; -$ztracker-100" | bc`
   zcalorimeter=`echo "scale=3; 160/2+30+1000" | bc`
#   divergence=`echo "scale=3; $i/100" | bc`
   emittance=`echo "scale=3; $i * 3 * 1.5 / 100" | bc`
   cpulimit -l 60 nice -n 19 Gate -a "'[energy,230] [dz,160] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,$i] [divergence,0.02] [emittance,$emittance] [rotation,0] [material,Water] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
   PIDLIST="$PIDLIST $!"
   echo "Running with beamspot = $i"

   if (( $IDX % $NCORES == 0)); then 
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

