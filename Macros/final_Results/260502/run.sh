#!/bin/bash
for sys in '0' '1'
do
root -l -b -q 'draw_Raa_JPsi_y0_1p6_pT_241014.C('$sys')'
root -l -b -q 'draw_Raa_JPsi_y1p6_2p4_pT_241014.C('$sys')'
root -l -b -q 'draw_Raa_JPsi_y0_1p6_Cent.C('$sys')'
root -l -b -q 'draw_Raa_JPsi_y1p6_2p4_Cent_4Bins_241014.C('$sys')'
root -l -b -q 'compare_pT_Jpsi.C('$sys')'
root -l -b -q 'compare_Npart_Jpsi.C('$sys')'
done
root -l -b -q 'compare_Bak_Gwak_Jpsi.C)'