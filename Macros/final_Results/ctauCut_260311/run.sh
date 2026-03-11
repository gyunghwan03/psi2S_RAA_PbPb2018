#!/bin/bash
root -l -b -q 'draw_Bfraction_y0_1p6.C'
root -l -b -q 'draw_Bfraction_y1p6_2p4.C'
root -l -b -q 'draw_Bfraction_y0_1p6_Jpsi.C'
root -l -b -q 'draw_Bfraction_y1p6_2p4_Jpsi.C'
for sys in '1'
do
	root -l -b -q 'draw_Raa_psi2S_y0_1p6_pT.C('$sys')' &
	root -l -b -q 'draw_Raa_psi2S_y1p6_2p4_pT.C('$sys')' &
	root -l -b -q 'draw_Raa_psi2S_y0_1p6_Cent.C('$sys')' &
	root -l -b -q 'draw_Raa_psi2S_y1p6_2p4_Cent_4Bins.C('$sys')' &

	root -l -b -q 'draw_Raa_JPsi_y0_1p6_Cent.C('$sys')' &
	root -l -b -q 'draw_Raa_JPsi_y0_1p6_pT_241014.C('$sys')' &
	root -l -b -q 'draw_Raa_JPsi_y1p6_2p4_Cent_4Bins_241014.C('$sys')' &
	root -l -b -q 'draw_Raa_JPsi_y1p6_2p4_pT_241014.C('$sys')' &
	root -l -b -q 'compare_Npart_Jpsi.C('$sys')' 
	root -l -b -q 'compare_pT_Jpsi.C('$sys')'
	root -l -b -q 'compare_pT_PR.C('$sys')' 
	root -l -b -q 'compare_pT_NP.C('$sys')' 
	root -l -b -q 'compare_Npart.C('$sys')' 
done

root -l -b -q 'DoubleRatio_Charmonia_Npart_midRap_240807.C' &
root -l -b -q 'DoubleRatio_Charmonia_Npart_fwdRap_240807.C' &
root -l -b -q 'DoubleRatio_Charmonia_pT_midRap.C' &
root -l -b -q 'DoubleRatio_Charmonia_pT_fwdRap.C'

