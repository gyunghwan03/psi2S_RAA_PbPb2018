#!/bin/bash

root -l -b <<EOF
.L Final2DFit.C 
.q
EOF

for pt in  '3,6.5' 
do
	root -l -b -q 'Final2DFit.C('$pt',1.6,2.4,0,180)'
done
for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,50'
do
	root -l -b -q 'Final2DFit.C('$pt',0,2.4,0,180)'
done
