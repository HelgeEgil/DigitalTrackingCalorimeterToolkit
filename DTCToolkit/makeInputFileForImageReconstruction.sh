#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

for y in `seq -w 0 5 0`; do
   root -l 'Scripts/makeInputFileForImageReconstruction.C(100, 800, '$y')'
done
