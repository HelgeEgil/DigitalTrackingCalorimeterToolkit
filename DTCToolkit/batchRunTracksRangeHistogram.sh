#!/bin/bash

echo "Usage: ./batchRunTracksRangeHistogram.sh <degraderThicknessFrom> <degraderThicknessStep> <degraderThicknessTo> <eventsPerRun> <runs>"

if [ $# -ne 5 ]; then
   echo "Incorrect number of arguments supplied, exiting"
   exit
fi

if [ $# -eq 5 ]; then
   IDX=0
   root -l -b -q 'Load.C'
   for i in `seq $1 $2 $3`;
   do
      IDX=$(( IDX+1 ))
      echo "Starting the $i mm degrader run with idx $IDX and absorber thickness $k."
      # Runs, dataType, recreate, energy, degraderthickness, eventsPerRun, outputFileIdx, drawFitResults, doTracking, excludeNuclearReactions, skipTracks
      tsp root -l -b -q 'Scripts/drawTracksRangeHistogramScript.C('$5', 0, 1, 917, '$i', '$4', '$IDX', true, true, false, 0)' # > OutputFiles/BPFitOutput_$IDX.txt
   done
   
   # Be sure that all jobs have completed, please tell me if there is a better way .... (tsp -w waits for the last job added to be completed, but not for the last job finished)
   for i in `seq $1 $2 $3`; do
      tsp -w
   done

   for j in `seq 1 1 $IDX`;
   do
      cat OutputFiles/result_makebraggpeakfit_idx$j.csv >> OutputFiles/result_makebraggpeakfit.csv
      rm OutputFiles/result_makebraggpeakfit_idx$j.csv
   done
fi
