#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

xlim=104

for spotx in `seq -84 4 104`; do
   tsp root -l 'Scripts/makeInputFileForCNNScript.C('$spotx')'
done
