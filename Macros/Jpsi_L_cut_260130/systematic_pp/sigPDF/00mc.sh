#!/bin/bash

root -l -b <<EOF
.L mc_MassFit_CBGauss.C 
.q
EOF

for pt in '3.5,6.5' '6.5,9' '9,12' '12,40' '3.5,40'
do
	root -l -b -q 'mc_MassFit_CBGauss.C('$pt',1.6,2.4)'
done
for pt in '6.5,9' '9,12' '12,15' '15,20' '20,25' '25,40' '6.5,40'
do
	root -l -b -q 'mc_MassFit_CBGauss.C('$pt',0,1.6)'
done

