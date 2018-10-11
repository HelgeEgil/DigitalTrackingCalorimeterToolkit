#!/bin/bash

NCORES=4
IDX=1

#for i in `seq -w 10 10 330`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
##   emittance=`echo "scale=3; 3*3 * 3" | bc`
#	Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,Water] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
#   PIDLIST="$PIDLIST $!"
#   echo "Running with i = $i"
#
#   if (( $IDX % $NCORES == 0)); then 
#      echo "Waiting for (PIDS $PIDLIST)"
#      time wait $PIDLIST
#      unset PIDLIST
#      IDX=1
#   else
#      IDX=$(( IDX+1 ))
#   fi
#done

for i in `seq -w 10 10 290`; do
   ztracker=`echo "scale=3; $i/2+15" | bc`
   zbeam=`echo "scale=3; -$ztracker-100" | bc`
   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
#   divergence=`echo "scale=3; $i/100" | bc`
#   emittance=`echo "scale=3; $i/100 * 9" | bc`
   Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,3] [emittance,15] [rotation,0] [material,A150] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
   PIDLIST="$PIDLIST $!"
   echo "Running with i = $i"

   if (( $IDX % $NCORES == 0 || ${i#0} == 290 )); then 
      echo "Waiting for (PIDS $PIDLIST)"
      time wait $PIDLIST
      unset PIDLIST
      IDX=1
   else
      IDX=$(( IDX+1 ))
   fi
done
#
#for i in `seq -w 10 10 200`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
##   emittance=`echo "scale=3; $i/100 * 9" | bc`
#	nice -n 19 cpulimit -l 60 Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,3] [emittance,15] [rotation,0] [material,CorticalBone] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
#   PIDLIST="$PIDLIST $!"
#   echo "Running with i = $i"
#
#   if (( $IDX % $NCORES == 0 || ${i#0} == 200 )); then 
#      echo "Waiting for (PIDS $PIDLIST)"
#      time wait $PIDLIST
#      unset PIDLIST
#      IDX=1
#   else
#      IDX=$(( IDX+1 ))
#   fi
#done
#
#for i in `seq -w 10 10 250`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
##   emittance=`echo "scale=3; $i/100 * 9" | bc`
#	nice -n 19 cpulimit -l 60 Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,3] [emittance,15] [rotation,0] [material,B100] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
#   PIDLIST="$PIDLIST $!"
#   echo "Running with i = $i"
#
#   if (( $IDX % $NCORES == 0 || ${i#0} == 250 )); then 
#      echo "Waiting for (PIDS $PIDLIST)"
#      time wait $PIDLIST
#      unset PIDLIST
#      IDX=1
#   else
#      IDX=$(( IDX+1 ))
#   fi
#done
#
#for i in `seq -w 10 10 350`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
##   emittance=`echo "scale=3; $i/100 * 9" | bc`
#	nice -n 19 cpulimit -l 60 Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,3] [emittance,15] [rotation,0] [material,myAdipose] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
#   PIDLIST="$PIDLIST $!"
#   echo "Running with i = $i"
#
#   if (( $IDX % $NCORES == 0 || ${i#0} == 350 )); then 
#      echo "Waiting for (PIDS $PIDLIST)"
#      time wait $PIDLIST
#      unset PIDLIST
#      IDX=1
#   else
#      IDX=$(( IDX+1 ))
#   fi
#done
