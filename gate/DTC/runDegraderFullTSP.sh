#!/bin/bash

NCORES=8
IDX=1

for i in `seq -w 160 1 160`; do
   Gate -a "'[energy,917] [degraderthickness,$i]'" Main_He_chip.mac
   PIDLIST="$PIDLIST $!"
   echo "Running: degrader thickness = $i mm"
done
