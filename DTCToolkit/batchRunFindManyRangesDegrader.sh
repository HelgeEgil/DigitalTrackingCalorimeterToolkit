#!/bin/bash

echo "Usage: ./batchRunTracksRangeHistogram.sh <degraderThicknessFrom> <degraderThicknessStep> <degraderThicknessTo> <absorberthicknessFrom> <absorberthicknessStep> <absorberthicknessTo>"

if [ $# -ne 6 ]; then
   echo "Incorrect number of arguments supplied, exiting"
   exit
fi

if [ $# -eq 6 ]; then
   for k in `seq $4 $5 $6`; do
      IDX=0
      for i in `seq $1 $2 $3`;
      do
         IDX=$(( IDX+1 ))
         echo "Starting the $i mm degrader run with idx $IDX and absorber thickness $k."
         if [ $IDX -eq 1 ]; then
            root -l -b -q 'Scripts/findManyRangesDegraderParallel.C('$k', '$i', '$IDX')'
            continue
         fi

         # Runs, dataType, recreate, energy, degraderthickness, eventsPerRun, outputFileIdx, drawFitResults, doTracking, excludeNuclearReactions, skipTracks
#         tsp root -l -b -q 'Scripts/findManyRangesDegraderParallel.C('$k', '$i', '$IDX')' 
      done
      
      # Be sure that all jobs have completed, please tell me if there is a better way .... (tsp -w waits for the last job added to be completed, but not for the last job finished)
#      for i in `seq $1 $2 $3`; do
#         tsp -w
#      done

      for j in `seq 1 1 $IDX`;
      do
         cat OutputFiles/findManyRangesDegrader_final_Helium_idx$j.csv >> OutputFiles/findManyRangesDegrader_final_Helium.csv
         rm OutputFiles/findManyRangesDegrader_final_Helium_idx$j.csv
      done
   done
fi
