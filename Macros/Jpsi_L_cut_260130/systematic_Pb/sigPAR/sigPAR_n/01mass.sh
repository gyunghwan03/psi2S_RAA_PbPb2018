#!/bin/bash

root -l -b <<EOF
.L MassFit_FixPar_Data.C 
.q
EOF

for PR in 1 2
do
	for pt in '3.5,6.5' '6.5,9' '9,12' '12,40'
	do
		root -l -b -q 'MassFit_FixPar_Data.C('$pt',1.6,2.4,0,180,'$PR')'
	done
	for pt in '6.5,9' '9,12' '12,15' '15,20' '20,25' '25,40'
	do
		root -l -b -q 'MassFit_FixPar_Data.C('$pt',0,1.6,0,180,'$PR')'
	done

	for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,180'
	do
		root -l -b -q 'MassFit_FixPar_Data.C(6.5,40,0,1.6,'$cent','$PR')'
	done
	for cent in '0,20' '20,60' '60,100' '100,180'
	do
		root -l -b -q 'MassFit_FixPar_Data.C(3.5,40,1.6,2.4,'$cent','$PR')'
	done
done