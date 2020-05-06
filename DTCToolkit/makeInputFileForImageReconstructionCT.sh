#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

for rot in `seq 0 2 358`; do
   if [ $rot != 90 ]; then
   tsp root -l 'Scripts/makeInputFileForImageReconstructionCT.C(250, 50, '$rot')'
   fi
done
