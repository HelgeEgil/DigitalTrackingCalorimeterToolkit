#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

for rot in `seq 0 2 358`; do
   tsp root -l 'Scripts/makeInputFileForImageReconstruction.C(1000, 100, '$rot')'
done
