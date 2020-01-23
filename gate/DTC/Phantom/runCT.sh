#!/bin/bash

rot="008"
echo "Running at rotation $rot deg"
for x in `seq -f "%04.f" -98 7 98`; do
   for y in `seq -f "%04.f" -70 7 70`; do
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
     tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[spotx,$x] [spoty,$y] [theta,$theta] [axisx,$axis_x] [axisy,$axis_y] [rotation,$rot]'" Main_chip_phantom.mac
   done
done

# Wait for all to finish here
tsp -w

# Combine all generated files into single file (1.5M protons)
rotNew=`echo "$rot"|bc`
root -l 'combineSpots.C('$rotNew')'

# Remove all spot-specific files
rm ../../../DTCToolkit/Data/MonteCarlo/DTC_Final_HeadPhantom_rotation${rot}deg_spot*

for rot in `seq -f "%03.0f" 12 2 358`; do
   echo "Running at rotation $rot deg"
   for x in `seq -f "%04.f" -98 7 98`; do
      for y in `seq -f "%04.f" -70 7 70`; do
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
        tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[spotx,$x] [spoty,$y] [theta,$theta] [axisx,$axis_x] [axisy,$axis_y] [rotation,$rot]'" Main_chip_phantom.mac
      done
   done

   # Wait for all to finish here
   tsp -w
   
   # Combine all generated files into single file (1.5M protons)
   rotNew=`echo "$rot"|bc`
   root -l 'combineSpots.C('$rotNew')'
 
   # Remove all spot-specific files
   rm ../../../DTCToolkit/Data/MonteCarlo/DTC_Final_HeadPhantom_rotation${rot}deg_spot*

done
