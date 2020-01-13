#!/bin/bash

for x in `seq -w -100 100 100`; do
   for y in `seq -w -100 10 100`; do
      arg=`echo "6600/sqrt($x^2+$y^2+6600^2)"|bc -l`
      theta=`echo "a(sqrt(1-1*$arg^2)/$arg)"|bc -l`
      if [ `echo "$x^2 + $y^2"|bc` -gt 0 ]; then
         axis_x=`echo "-1*$y/sqrt($x^2+$y^2)"|bc -l`
         axis_y=`echo "$x/sqrt($x^2+$y^2)"|bc -l`
      else
         axis_x=1
         axis_y=0
      fi
   #   tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[spotx,$x] [spoty,$y] [theta,$theta] [axisx,$axis_x] [axisy,$axis_y] [rotation,90]'" Main_chip_phantom.mac
      /home/rttn/gate/Gate-8.2custom-install/bin/Gate -a "'[spotx,$x] [spoty,$y] [theta,$theta] [axisx,$axis_x] [axisy,$axis_y] [rotation,90]'" Main_chip_phantom.mac
      PIDLIST="$PIDLIST $!"
      echo "Running at spot ($x,$y) -- total angle $theta and axis of rotation $axis_x $axis_y"
   done
done
