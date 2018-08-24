#!/bin/bash
# Usage: drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>

echo "Usage: ./drawTracks3D.sh <eventsPerRun> <runs> <degraderthickness>"

if [ $# -ne 3 ]; then
	echo "Incorrect number of arguments supplied, exiting"
   echo eventsPerRun: $1, runs: $2, degraderthickness: $3
	exit
fi

echo eventsPerRun: $1, runs: $2, degraderthickness: $3

root -l 'Scripts/drawTracks3DScript.C('$2', 0, 1, 100, 250, '$3', '$1', true)'
