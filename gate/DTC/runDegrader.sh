#!/bin/bash

NCORES=8
IDX=1

for i in `seq -w 160 1 160`; do
#   hdt=`echo "scale=3; -$i/2-2" | bc`
   Gate -a "'[energy,230] [degraderthickness,$i]'" Main.mac
   PIDLIST="$PIDLIST $!"
   echo "Running: degrader thickness = $i mm"
done
