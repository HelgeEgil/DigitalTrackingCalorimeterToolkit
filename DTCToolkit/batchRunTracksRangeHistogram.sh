#!/bin/bash

echo "Usage: ./batchRunTracksRangeHistogram.sh <degraderThicknessFrom> <degraderThicknessStep> <degraderThicknessTo> <eventsPerRun> <runs> <absorberthicknessFrom> <absorberthicknessStep> <absorberthicknessTo>"

if [ $# -ne 8 ]; then
   echo "Incorrect number of arguments supplied, exiting"
   exit
fi

if [ $# -eq 8 ]; then
   for k in `seq $6 $7 $8`; do
      IDX=0
      # Change Absorber Thickness
      echo "Changing absorber thickness to $k..."
      sed -i "s/kAbsorberThickness = [0-9];/kAbsorberThickness = $k;/" GlobalConstants/Constants.h
      # COMPILE
      root -l -b -q 'Load.C'
      for i in `seq $1 $2 $3`;
      do
         IDX=$(( IDX+1 ))
         echo "Starting the $i mm degrader run with idx $IDX and absorber thickness $k."
         # Runs, dataType, recreate, energy, degraderthickness, eventsPerRun, outputFileIdx, drawFitResults, doTracking, excludeNuclearReactions, skipTracks
         root -l -b -q 'Scripts/drawTracksRangeHistogramScript.C('$5', 0, 1, 230, '$i', '$4', '$IDX', false, false, false, 0)' > OutputFiles/BPFitOutput_$IDX.txt &
      done
      
      # Be sure that all jobs have completed, please tell me if there is a better way .... (tsp -w waits for the last job added to be completed, but not for the last job finished)
#      for i in `seq $1 $2 $3`; do
#         tsp -w
#      done

      for j in `seq 1 1 $IDX`;
      do
         cat OutputFiles/result_makebraggpeakfit_idx$j.csv >> OutputFiles/result_makebraggpeakfit.csv
         rm OutputFiles/result_makebraggpeakfit_idx$j.csv
      done
   done
fi
