#!/bin/bash

root -l -b <<EOF
.L MassFit_FixPar_Data.C 
.q
EOF

#for pt in  '3,4.5' '4.5,6.5' 
for pt in  '3,6.5'
do
	root -l -b -q 'MassFit_FixPar_Data_LowPt.C('$pt',1.6,2.4)'
done
#for pt in  '6.5,7' '7,7.5' '7.5,8' '8,9' '9,10' '10,12' '12,15' '15,20' '20,50' '6.5,9' '9,12' 
for pt in '6.5,9' '9,12' '12,15' '15,20' '20,50'
do
	root -l -b -q 'MassFit_FixPar_Data.C('$pt',0,1.6)'
done
#for pt in '6.5,8' '8,12' '12,50' '6.5,12'
for pt in '6.5,12' '12,50'
#for pt in '6.5,8' '8,9' '9,10' '10,12' '12,15'
do
	root -l -b -q 'MassFit_FixPar_Data.C('$pt',1.6,2.4)'
done

#for pt in '6.5,8' '8,9' '9,10' '10,11' '11,12' '12,13.5' '13.5,15'
#do
#	root -l -b -q 'MassFit_FixPar_Data.C('$pt',0,1.2)'
#done
