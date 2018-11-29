#!/bin/bash

NCORES=4
IDX=1

for i in `seq -f "%03g" 160 10 160`; do
   ztracker=`echo "scale=3; $i/2+15" | bc`
   zbeam=`echo "scale=3; -$ztracker-100" | bc`
   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
#   divergence=`echo "scale=3; $i/100" | bc`
   emittance=`echo "scale=3; $i * 3 * 1.5" | bc`
   cpulimit -l 60 nice -n 60 Gate -a "'[energy,917] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,Water] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
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

echo "Waiting for (PIDS $PIDLIST)"
time wait $PIDLIST
unset PIDLIST

#for i in `seq -f "%03g" 10 10 40`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
#   emittance=`echo "scale=3; $i * 3 * 1.5" | bc`
#   nice -n 19 cpulimit -l 60 Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,A150] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
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
#
#echo "Waiting for (PIDS $PIDLIST)"
#time wait $PIDLIST
#unset PIDLIST
#
#for i in `seq -f "%03g" 10 10 40`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
#   emittance=`echo "scale=3; $i * 3 * 1.5" | bc`
#	nice -n 19 cpulimit -l 60 Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,CorticalBone] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
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
#
#echo "Waiting for (PIDS $PIDLIST)"
#time wait $PIDLIST
#unset PIDLIST
#
#for i in `seq -f "%03g" 10 10 40`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
#   emittance=`echo "scale=3; $i * 3 * 1.5" | bc`
#	nice -n 19 cpulimit -l 60 Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,B100] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
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
#
#echo "Waiting for (PIDS $PIDLIST)"
#time wait $PIDLIST
#unset PIDLIST
#
#for i in `seq -f "%03g" 10 10 40`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zbeam=`echo "scale=3; -$ztracker-100" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
##   divergence=`echo "scale=3; $i/100" | bc`
#   emittance=`echo "scale=3; $i * 3 * 1.5" | bc`
#	nice -n 19 cpulimit -l 60 Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [zbeam,$zbeam] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,myAdipose] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
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
#
#echo "Waiting for (PIDS $PIDLIST)"
#time wait $PIDLIST
#unset PIDLIST
