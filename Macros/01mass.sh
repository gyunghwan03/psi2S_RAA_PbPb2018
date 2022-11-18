#!/bin/bash

root -l -b <<EOF
.L MassFit_FixPar_Data.C 
.q
EOF

for pt in  '3,6.5' '6.5,12' '12,50'
do
	root -l -b -q 'MassFit_FixPar_Data.C('$pt',1.6,2.4,0,200)'
done
for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,50'
do
	root -l -b -q 'MassFit_FixPar_Data.C('$pt',0,1.6,0,200)'
done
