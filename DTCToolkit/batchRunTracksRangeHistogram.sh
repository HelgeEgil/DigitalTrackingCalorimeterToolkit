#!/bin/bash

NCORES=8
IDX=1

echo "Usage: ./batchRunTracksRangeHistogram.sh <degraderThicknessFrom> <degraderThicknessStep> <degraderThicknessTo> <eventsPerRun> <runs>"

if [ $# -ne 5 ]; then
   echo "Incorrect number of arguments supplied, exiting"
   exit
fi

if [ $# -eq 5 ]; then
   # COMPILE BEFORE ANYTHING
   root -l -q 'Load.C'
   for i in `seq $1 $2 $3`;
   do
      echo "Making plots at $i mm degrader with idx $IDX."
      root -l -q 'Scripts/drawTracksRangeHistogramScript.C('$5', 0, 1, 250, '$i', '$4', '$IDX', false, true, false)' > OutputFiles/BPFitOutput_$IDX.txt &
      PIDLIST="$PIDLIST $!"

      if (( $IDX % $NCORES == 0 || $i == $3)); then
         echo "Waiting for all started runs... (PIDS $PIDLIST)"
         time wait $PIDLIST
         unset PIDLIST

         echo "Concatenating result_makebraggpeakfit files before next loop"
         for j in `seq 1 1 $IDX`;
         do
            cat OutputFiles/result_makebraggpeakfit_idx$j.csv >> OutputFiles/result_makebraggpeakfit.csv
            rm OutputFiles/result_makebraggpeakfit_idx$j.csv
         done
         IDX=1

      else
         echo "Increasing IDX"
         IDX=$(( IDX+1 ))
      fi
   done
fi
