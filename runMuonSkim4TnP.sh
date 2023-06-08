#!/bin/bash


# Mid Rapidity Prompt
	# Nominal
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 0)'

	# ID Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 1, 0, 0)'
	# ID Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 2, 0, 0)'
	# ID Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, -1, 0, 0)'
	# ID Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, -2, 0, 0)'
	# ID Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 99, 0, 0)'

	# Trk Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 1, 0)'
	# Trk Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 2, 0)'
	# Trk Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, -1, 0)'
	# Trk Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, -2, 0)'
	# Trk Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 99, 0)'

	# Trg Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 1)'
	# Trg Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 2)'
	# Trg Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -1)'
	# Trg Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -2)'
	# Trg Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 99)'

# Forward Rapidity Prompt
	# Nominal
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 0)'

	# ID Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 1, 0, 0)'
	# ID Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 2, 0, 0)'
	# ID Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, -1, 0, 0)'
	# ID Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, -2, 0, 0)'
	# ID Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 99, 0, 0)'

	# Trk Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 1, 0)'
	# Trk Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 2, 0)'
	# Trk Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, -1, 0)'
	# Trk Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, -2, 0)'
	# Trk Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 99, 0)'

	# Trg Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 1)'
	# Trg Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 2)'
	# Trg Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -1)'
	# Trg Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -2)'
	# Trg Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, true, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 99)'





# Mid Rapidity Non-prompt
	# Nominal
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 0)'

exit

	# ID Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 1, 0, 0)'
	# ID Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 2, 0, 0)'
	# ID Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, -1, 0, 0)'
	# ID Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, -2, 0, 0)'
	# ID Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 99, 0, 0)'

	# Trk Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 1, 0)'
	# Trk Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 2, 0)'
	# Trk Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, -1, 0)'
	# Trk Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, -2, 0)'
	# Trk Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 99, 0)'

	# Trg Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 1)'
	# Trg Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 2)'
	# Trg Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -1)'
	# Trg Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -2)'
	# Trg Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 0.0, 1.6, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 99)'

# Forward Rapidity Non-prompt
	# Nominal
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 0)'

	# ID Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 1, 0, 0)'
	# ID Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 2, 0, 0)'
	# ID Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, -1, 0, 0)'
	# ID Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, -2, 0, 0)'
	# ID Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 99, 0, 0)'

	# Trk Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 1, 0)'
	# Trk Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 2, 0)'
	# Trk Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, -1, 0)'
	# Trk Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, -2, 0)'
	# Trk Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 99, 0)'

	# Trg Stat Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 1)'
	# Trg Stat Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 2)'
	# Trg Syst Up
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -1)'
	# Trg Syst Do
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, -2)'
	# Trg Tag Change
root -l -b -q 'makeMuonSkimTree.C(true, false, 0.0, 50.0, 1.6, 2.4, 0, 180, true, true, false, kTrigJpsi, 0, 0, 0, 99)'
