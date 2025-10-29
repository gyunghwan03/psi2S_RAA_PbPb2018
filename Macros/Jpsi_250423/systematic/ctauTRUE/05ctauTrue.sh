#!/bin/bash

root -l -b <<EOF
.L CtauTrue.C 
.q
EOF

for pt in  '3.5,6.5' '6.5,9' '9,12' '12,40'
do
	root -l -b -q 'CtauTrue.C('$pt',1.6,2.4,0,180)'
done
for pt in  '6.5,9' '9,12' '12,15' '15,20' '20,25' '25,40'
do
	root -l -b -q 'CtauTrue.C('$pt',0,1.6,0,180)'
done

for cent in '0,180'
do
    root -l -b -q 'CtauTrue.C(6.5,40,0,1.6,'$cent')'
done
for cent in '0,180'
do
    root -l -b -q 'CtauTrue.C(3.5,40,1.6,2.4,'$cent')'
done
