#!/bin/sh

make || exit

if [ $# -eq 1 ]
then 
    sh testDijkstra.sh $1
    sh testCentralityMeasures.sh $1
    sh testLanceWilliamsHAC.sh $1
else
    sh testDijkstra.sh 
    sh testCentralityMeasures.sh
    sh testLanceWilliamsHAC.sh
fi
exit 1