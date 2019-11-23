#!/bin/bash

NCORES=8
IDX=1

for i in `seq -w 274 4 330`; do
   hdt=`echo "scale=3; -$i/2-2" | bc`
   beampos=`echo "scale=3; -$i-5" | bc`
   tsp Gate -a "'[energy,917] [degraderthickness,$i] [halfdegraderthickness,$hdt] [beampos,$beampos]'" Main.mac
   PIDLIST="$PIDLIST $!"
   echo "Running: Absorber thickness = $j mm; degrader thickness = $i mm"
done
