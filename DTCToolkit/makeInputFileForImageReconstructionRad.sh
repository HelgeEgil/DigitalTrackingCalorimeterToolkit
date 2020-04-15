#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

for spotx in `seq -98 7 98`; do
   tsp root -l 'Scripts/makeInputFileForImageReconstructionRad.C(500, 50, 90, '$spotx')'
done
