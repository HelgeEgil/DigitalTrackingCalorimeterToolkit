#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

if [ $# -ne 1 ]; then
   echo "Please give phantom as argument (linePair, CTP404, headphantom, wedge)"
   exit
fi

xlim=104

if [ $1 == "headphantom" ]; then
   xlim=98
fi

for spotx in `seq -$xlim 7 $xlim`; do
#   tsp root -l 'Scripts/makeInputFileForImageReconstructionRad.C(500, 50, 90, '$spotx', "'$1'")'
   tsp root -l 'Scripts/makeInputFileForImageReconstructionRad.C(200, 50, 90, '$spotx', "'$1'")'
done
