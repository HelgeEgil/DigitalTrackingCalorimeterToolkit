#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

for rot in `seq 0 4 358`; do
   tsp root -l 'Scripts/makeInputFileForImageReconstructionCT.C(150, 100, '$rot')'
done
