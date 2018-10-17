#!/bin/bash

NCORES=4
IDX=1

for i in `seq -w 200 1 200`; do
   ztracker=`echo "scale=3; $i/2+15" | bc`
   zbeam=`echo "scale=3; -$ztracker-100" | bc`
   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
#   divergence=`echo "scale=3; $i/100" | bc`
#   emittance=`echo "scale=3; $i * 3 * 1.5" | bc`
	nice -n 19 cpulimit -l 60 Gate -a "'[energy,200] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,Water] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
   PIDLIST="$PIDLIST $!"
   echo "Running with i = $i"

   if (( $IDX % $NCORES == 0)); then 
      echo "Waiting for (PIDS $PIDLIST)"
      time wait $PIDLIST
      unset PIDLIST
      IDX=1
   else
      IDX=$(( IDX+1 ))
   fi
done
