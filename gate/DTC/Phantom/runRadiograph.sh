#!/bin/bash

#phantom="CTP404" # x = `seq -77 7 77`; y = `seq -21 7 21`
#phantom="linePair" # x =`seq -77 7 77`; y = `seq -21 7 21`
#phantom="headphantom" # x = `seq -98 7 98`; y = `seq -88 7 88`
phantom="wedge"

step=4
xlim=84
ylim=28

if [ $phantom == "headphantom" ]; then
   xlim=98
   ylim=88
fi

if [ $phantom == "wedge" ]; then
   xlim=104
   ylim=0
fi


for x in `seq -f "%04.f" -84 $step $xlim`; do
   for y in `seq -f "%04.f" 0 $step $ylim`; do
      arg=`echo "6600/sqrt($x^2+$y^2+6600^2)"|bc -l`
      theta=`echo "a(sqrt(1-1*$arg^2)/$arg)"|bc -l`
      if [ `echo "$x^2 + $y^2"|bc` -gt 0 ]; then
         axis_x=`echo "-1*$y/sqrt($x^2+$y^2)"|bc -l`
         axis_y=`echo "$x/sqrt($x^2+$y^2)"|bc -l`
      else
         axis_x=1
         axis_y=0
      fi
     tsp /build/gate/Gate-8.2reducedRootOutput-install/bin/Gate -a "'[phantom,$phantom] [spotx,$x] [spoty,$y] [theta,$theta] [axisx,$axis_x] [axisy,$axis_y] [rotation,90]'" Main_phantom.mac
      PIDLIST="$PIDLIST $!"
      echo "Running at spot ($x,$y) -- total angle $theta and axis of rotation $axis_x $axis_y"
   done
done
