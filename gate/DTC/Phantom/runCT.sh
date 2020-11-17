#!/bin/bash

#phantom="CTP404" # x = `seq -77 7 77`; y = `seq -21 7 21`
# phantom="linePair" # x =`seq -77 7 77`; y = `seq -21 7 21`
#phantom="headphantom" # x = `seq -98 7 98`; y = `seq -88 7 88`
#
#for rot in `seq -f "%03.0f" 0 2 358`; do
#   echo "Running at rotation $rot deg"
#   for x in `seq -f "%04.f" -98 7 98`; do
#      for y in `seq -f "%04.f" -88 7 88`; do
#         arg=`echo "6600/sqrt($x^2+$y^2+6600^2)"|bc -l`
#         theta=`echo "a(sqrt(1-1*$arg^2)/$arg)"|bc -l`
#         if [ `echo "$x^2 + $y^2"|bc` -gt 0 ]; then
#            axis_x=`echo "-1*$y/sqrt($x^2+$y^2)"|bc -l`
#            axis_y=`echo "$x/sqrt($x^2+$y^2)"|bc -l`
#         else
#            # Don't calculate axis of rotation if theta==0
#            axis_x=1
#            axis_y=0
#         fi
#        tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[phantom,$phantom] [rotation,$rot] [spotx,$x] [spoty,$y] [theta,$theta] [axisx,$axis_x] [axisy,$axis_y]'" Main_phantom.mac
#      done
#   done
#
#   # Wait for all to finish here
#   tsp -w
#   
#   # Combine all generated files into single file (1.5M protons)
#   rotNew=`echo "$rot"|bc`
#   root -l 'combineSpots_'${phantom}'.C('$rotNew')'
# 
#   # Remove all spot-specific files
#   rm ../../../DTCToolkit/Data/MonteCarlo/DTC_Final_${phantom}_rotation${rot}deg_spot*
#done

phantom="CTP404" # x = `seq -77 7 77`; y = `seq -21 7 21`
#phantom="linePair" # x =`seq -77 7 77`; y = `seq -21 7 21`
# phantom="headphantom" # x = `seq -98 7 98`; y = `seq -88 7 88`
#
for rot in `seq -f "%03.0f" 180 1 359`; do
   echo "Running at rotation $rot deg"
   for x in `seq -f "%04.f" -84 7 84`; do
      for y in `seq -f "%04.f" -28 7 28`; do
         arg=`echo "6600/sqrt($x^2+$y^2+6600^2)"|bc -l`
         theta=`echo "a(sqrt(1-1*$arg^2)/$arg)"|bc -l`
         if [ `echo "$x^2 + $y^2"|bc` -gt 0 ]; then
            axis_x=`echo "-1*$y/sqrt($x^2+$y^2)"|bc -l`
            axis_y=`echo "$x/sqrt($x^2+$y^2)"|bc -l`
         else
            # Don't calculate axis of rotation if theta==0
            axis_x=1
            axis_y=0
         fi
        tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[phantom,$phantom] [rotation,$rot] [spotx,$x] [spoty,$y] [theta,$theta] [axisx,$axis_x] [axisy,$axis_y]'" Main_phantom.mac
      done
   done

   # Wait for all to finish here
   tsp -w
   
   # Combine all generated files into single file (1.5M protons)
   rotNew=`echo "$rot"|bc`
   root -l 'combineSpots_'${phantom}'.C('$rotNew')'
 
   # Remove all spot-specific files
   rm ../../../DTCToolkit/Data/MonteCarlo/DTC_Final_${phantom}_rotation${rot}deg_spot*
done
