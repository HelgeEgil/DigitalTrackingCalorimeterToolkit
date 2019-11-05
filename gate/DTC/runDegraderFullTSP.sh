#!/bin/bash

NCORES=8
IDX=1

for j in `seq -w 35 1 35`; do
   for i in `seq -w 5 5 330`; do
      hdt=`echo "scale=3; -$i/2-2" | bc`
      beampos=`echo "scale=3; -$i-5" | bc`
      cpulimit -l 70 Gate -a "'[absorberthickness,$j] [energy,917] [degraderthickness,$i] [halfdegraderthickness,$hdt] [beampos,$beampos]'" Main.mac > terminal_output.txt &
      PIDLIST="$PIDLIST $!"
      echo "Running (full): Absorber thickness = $j mm; degrader thickness = $i mm"
   done
done
