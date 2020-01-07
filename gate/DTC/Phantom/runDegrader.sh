#!/bin/bash

for i in `seq -w 2 4 330`; do
#   hdt=`echo "scale=3; -$i/2-2" | bc`
   tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[energy,230] [degraderthickness,$i]'" Main_He_full.mac
   PIDLIST="$PIDLIST $!"
   echo "Running: degrader thickness = $i mm"
done
