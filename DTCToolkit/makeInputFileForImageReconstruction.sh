#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

for y in `seq -w -75 5 75`; do
   tsp root -l 'Scripts/makeInputFileForImageReconstruction.C(100, 800, $y)'
done
