#!/bin/bash
# Usage: run.sh material energy_from energy_increment energy_to

NCORES=3
IDX=1

echo "Usage: ./makeBPFitResults.sh <degraderthickness_from> <degraderthickness_increment> <degrader_thickness_to> <absorberthickness_from> <absorberthickness_increment> <absorberthickness_to>"
echo degraderthickness_from: $1, degraderthickness_increment: $2 degraderthickness_to: $3

if [ $# -eq 6 ]; then
	echo "Looping over absorberthicknesses!"
   echo absorberthickness_from: $4, absorberthickness_increment: $5 absorberthickness_to: $6

   for j in `seq $4 $5 $6`;
   do
      sed -i "58s/[2-6]/$j/" GlobalConstants/Constants.h
      # COMPILE BEFORE ANYTHING
      root -l -q 'Load.C'

      for i in `seq $1 $2 $3`;
      do
      	echo "Making plots at $i mm degrader with $j mm absorber and idx $IDX."
	      root -l -q 'Scripts/makeBraggPeakPDFDegrader.C('$i,$IDX')' > OutputFiles/BPFitOutput_$IDX.txt &
         PIDLIST="$PIDLIST $!"
      
         if (( $IDX % $NCORES == 0 || $i == $3)); then
            echo "Waiting for all started runs... (PIDS $PIDLIST)"
            time wait $PIDLIST
            unset PIDLIST

            echo "Concatenating result_makebraggpeakfit files before next loop"
            for k in `seq 1 1 $IDX`;
            do
               cat OutputFiles/result_makebraggpeakfit_idx$k.csv >> OutputFiles/result_makebraggpeakfit.csv
               rm OutputFiles/result_makebraggpeakfit_idx$k.csv
            done
            IDX=1

         else
            echo "Increasing IDX"
            IDX=$(( IDX+1 ))
         fi
      done
   done
fi


if [ $# -eq 3 ]; then
   echo "Running with absorberthickness in Constants.h"
   # COMPILE BEFORE ANYTHING
   root -l -q 'Load.C'
   for i in `seq $1 $2 $3`;
   do
      echo "Making plots at $i mm degrader with idx $IDX."
      root -l -q 'Scripts/makeBraggPeakPDFDegrader.C('$i,$IDX')' > OutputFiles/BPFitOutput_$IDX.txt &
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
