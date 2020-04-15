#!/bin/bash

for i in `seq -w 1 50 301`; do
   Gate -a "'[energy,230] [degraderthickness,$i]'" Main_chip.mac
   echo "Running PROTON degrader thickness = $i mm"
done

#for i in `seq -w 50 1 50`; do
   # tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[energy,917] [degraderthickness,$i]'" Main_He_chip.mac
#   ~/gate/Gate-8.2custom-install/bin/Gate -a "'[energy,917] [degraderthickness,$i]'" Main_He_chip.mac
#   Gate -a "'[energy,917] [degraderthickness,$i]'" Main_He_chip.mac
#   echo "Running HELIUM degrader thickness = $i mm"
#done
