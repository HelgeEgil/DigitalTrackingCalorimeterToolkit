#!/bin/bash
# Usage: drawTracksDepthDose.sh <eventsPerRun> <runs> <degraderthickness>

echo "Usage: ./drawTracksDepthDose.sh <eventsPerRun> <runs> <degraderthickness> "

if [ $# -ne 3 ]; then
	echo "Incorrect number of arguments supplied, exiting"
	exit
fi

echo eventsPerRun: $1, runs: $2, degraderthickness: $3

# Arguments: Runs, Use MC data, Remake tracks, Initial energy, Phantom degrader thickness (file to use), Events Per Run, File numbering for batch jobs (-1 = default), draw histogram, perform tracking (instead of using MC truth), removed MC-tagged tracks with nuclear events
root -l 'Scripts/drawTracksRangeHistogramScript.C('$2', 0, 1, 250, '$3', '$1', -1, true, true, false)'
