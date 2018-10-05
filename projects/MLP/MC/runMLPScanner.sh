#!/bin/bash

echo "Usage: ./run.sh <energy> <phantomthickness_from> <phantomthickness_increment> <phantomthickness_to>"
echo Initial Energy: $1 [MeV], Phantom thickness [mm] from $2 step $3 to $4

NCORES=4
IDX=1

#if [ $# -ne 4 ]; then
#	echo Invalid number of arguments: $#
#	exit
#fi

for i in `seq 10 10 320`; do
   ztracker=`echo "scale=3; $i/2+15" | bc`
   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
#   divergence=`echo "scale=3; $i/100" | bc`
#   emittance=`echo "scale=3; $i/100 * 9" | bc`
	Gate -a "'[energy,230] [dz,$i] [ztracker,$ztracker] [spotsize,3] [divergence,2] [emittance,15] [rotation,0] [material,Water] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
  # Gate -a "'[energy,230] [dz,160] [spotsize,3] [divergence,$divergence] [emittance,$emittance] [ztracker,$ztracker] [zcalorimeter,$zcalorimeter] [material,]'" MLPscanner.mac  > terminal_output.txt &
   PIDLIST="$PIDLIST $!"
   echo "Running with i = $i"

   if (( $IDX % $NCORES == 0 || $i == 50 )); then 
      echo "Waiting for (PIDS $PIDLIST)"
      time wait $PIDLIST
      unset PIDLIST
      IDX=1
   else
      IDX=$(( IDX+1 ))
   fi
done

#IDX=1
#for i in `seq 20 20 160`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
#	#nice -n 19 cpulimit -l 50 Gate -a "'[energy,$1] [dz,$i] [ztracker,$ztracker] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
#   nice -n 19 cpulimit -l 50 Gate -a "'[energy,200] [dz,$i] [spotsize,3] [ztracker,$ztracker] [zcalorimeter,$zcalorimeter] [material,CorticalBone]'" MLPscanner.mac &
#   PIDLIST="$PIDLIST $!"
#   echo "Running with spotsize $spotsize mm"
#
#   if (( $IDX % $NCORES == 0 || $i == 160 )); then 
#      echo "Waiting for (PIDS $PIDLIST)"
#      time wait $PIDLIST
#      unset PIDLIST
#      IDX=1
#   else
#      IDX=$(( IDX+1 ))
#   fi
#done
#
#IDX=1
#for i in `seq 20 20 230`; do
#   ztracker=`echo "scale=3; $i/2+15" | bc`
#   zcalorimeter=`echo "scale=3; $i/2+30+1000" | bc`
#	#nice -n 19 cpulimit -l 50 Gate -a "'[energy,$1] [dz,$i] [ztracker,$ztracker] [zcalorimeter,$zcalorimeter]'" MLPscanner.mac > terminal_output.txt &
#   nice -n 19 cpulimit -l 50 Gate -a "'[energy,200] [dz,$i] [spotsize,3] [ztracker,$ztracker] [zcalorimeter,$zcalorimeter] [material,A150]'" MLPscanner.mac &
#   PIDLIST="$PIDLIST $!"
#   echo "Running with spotsize $spotsize mm"
#
#   if (( $IDX % $NCORES == 0 || $i == 230 )); then 
#      echo "Waiting for (PIDS $PIDLIST)"
#      time wait $PIDLIST
#      unset PIDLIST
#      IDX=1
#   else
#      IDX=$(( IDX+1 ))
#   fi
#done
#
