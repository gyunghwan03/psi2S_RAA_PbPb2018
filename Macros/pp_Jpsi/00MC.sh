#!/bin/bash

root -l -b <<EOF
.L mc_MassFit_HighpT.C 
.q
EOF

#for pt in  '3.5,4.5' '4.5,6.5' '6.5,7' '7,8' '8,10' '10,12' '12,15' '15,20' '20,50'
for pt in  '3,6.5' '6.5,9' '9,12' '12,50'
do
	root -l -b -q 'mc_MassFit_HighpT.C('$pt',1.6,2.4)'
done
#for pt in  '6.5,7' '7,7.5' '7.5,8' '8,9' '9,10' '10,12' '12,15' '15,20' '20,50' '6.5,9' '9,12' 
for pt in '6.5,9' '9,12' '12,15' '15,20' '20,25' '25,50'
do
	root -l -b -q 'mc_MassFit_HighpT.C('$pt',0,1.6)'
done
