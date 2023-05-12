#!/bin/bash

root -l -b <<EOF
.L Final2DFit.C 
.q
EOF

for y in '0,0.4' '0.4,0.8' '0.8,1.2' '1.2,1.6' '1.6,2.0' '2.0,2.4'
do
    root -l -b -q 'Final2DFit.C(6.5,50,'$y')'
done

#for pt in  '3,6.5' '6.5,12' '12,50' 
#for pt in '3,4' '4,5' '5,6' 
#for pt in '3,6.5' 
#do
#	root -l -b -q 'Final2DFit_LowPt.C('$pt',1.6,2.4)'
#done
#for pt in '6.5,12' '12,50' 
#do
#	root -l -b -q 'Final2DFit.C('$pt',1.6,2.4)'
#done
#for pt in '6.5,9' '9,12' '12,15' '15,20' '20,50'
#do
#	root -l -b -q 'Final2DFit.C('$pt',0,1.6)'
#done
#for pt in '6.5,8' '8,9' '9,10' '10,12' '12,15'
#do
#	root -l -b -q 'Final2DFit.C('$pt',1.6,2.4)'
#done
#for pt in '6.5,8' '8,9' '9,10' '10,11' '11,12' '12,13.5' '13.5,15'
#do
#	root -l -b -q 'Final2DFit.C('$pt',0,1.2)'
#done
