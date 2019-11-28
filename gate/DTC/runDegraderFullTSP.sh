#!/bin/bash

NCORES=8
IDX=1

for i in `seq -w 160 1 160`; do
   hdt=`echo "scale=3; -$i/2-2" | bc`
   beampos=`echo "scale=3; -$i-5" | bc`
   tsp Gate -a "'[energy,917] [degraderthickness,$i] [halfdegraderthickness,$hdt] [beampos,$beampos]'" Main_He.mac
   PIDLIST="$PIDLIST $!"
   echo "Running: degrader thickness = $i mm"
done
