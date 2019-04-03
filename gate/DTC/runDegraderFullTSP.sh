#!/bin/bash

NCORES=4
IDX=1

for j in `seq -w 2 1 6`; do
   for i in `seq -w 1 1 370`; do
      hdt=`echo "scale=3; -$i/2-2" | bc`
      beampos=`echo "scale=3; -$i-5" | bc`
      tsp Gate -a "'[absorberthickness,$j] [energy,250] [degraderthickness,$i] [halfdegraderthickness,$hdt] [beampos,$beampos]'" Main_full.mac > terminal_output.txt &
      PIDLIST="$PIDLIST $!"
      echo "Running (full): Absorber thickness = $j mm; degrader thickness = $i mm"
   done
done
