#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

for x in `seq -100 5 100`; do
   tsp root -l 'Scripts/makeInputFileForImageReconstruction.C(150, 100, '$x')'
done
