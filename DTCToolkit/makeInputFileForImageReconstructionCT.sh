#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

if [ $# -ne 1 ]; then
   echo "Please give phantom as argument (linePair, CTP404, headphantom, wedge)"
   exit
fi

for rot in `seq 0 1 359`; do
   tsp root -l 'Scripts/makeInputFileForImageReconstructionCT.C(200, 50, '$rot', "'$1'")'
done
