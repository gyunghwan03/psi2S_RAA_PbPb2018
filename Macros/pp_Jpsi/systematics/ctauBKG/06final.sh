#!/bin/bash

root -l -b <<EOF
.L Final2DFit.C 
.q
EOF

#for y in '0,0.4' '0.4,0.8' '0.8,1.2' '1.2,1.6' '1.6,2.0' '2.0,2.4'
#do
#    root -l -b -q 'Final2DFit.C(6.5,50,'$y')'
#done

for pt in  '3.5,6.5' '6.5,9' '9,12' '12,40' '3.5,40' 
do
	root -l -b -q 'Final2DFit.C('$pt',1.6,2.4)'
done
for pt in '6.5,9' '9,12' '12,15' '15,20' '20,25' '25,40' '6.5,40'
do
	root -l -b -q 'Final2DFit.C('$pt',0,1.6)'
done
#for pt in '6.5,8' '8,9' '9,10' '10,12' '12,15'
#do
#	root -l -b -q 'Final2DFit.C('$pt',1.6,2.4)'
#done
#for pt in '6.5,8' '8,9' '9,10' '10,11' '11,12' '12,13.5' '13.5,15'
#do
#	root -l -b -q 'Final2DFit.C('$pt',0,1.2)'
#done
