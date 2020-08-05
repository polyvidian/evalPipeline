import sys
from rdkit import Chem
from rdkit.Chem import rdinchi, Descriptors


canDeriv = set()
tmsTargetInchis = set()
canDerivInchis = set()

theOH = Chem.MolFromSmarts('[O,o][H]')
theSH = Chem.MolFromSmarts('[S,s][H]')
theNH = Chem.MolFromSmarts('[N,n][H]')
theTMS = Chem.MolFromSmiles("[Si](C)(C)C")

Cpatt = Chem.MolFromSmarts("[C,c]")


def FragIndicesToMol(oMol, indices):
	em = Chem.EditableMol(Chem.Mol())

	newIndices = {}
	for i, idx in enumerate(indices):
		em.AddAtom(oMol.GetAtomWithIdx(idx))
		newIndices[idx] = i

	for i, idx in enumerate(indices):
		at = oMol.GetAtomWithIdx(idx)
		for bond in at.GetBonds():
			if bond.GetBeginAtomIdx() == idx:
				oidx = bond.GetEndAtomIdx()
			else:
				oidx = bond.GetBeginAtomIdx()
			# make sure every bond only gets added once:
			if oidx < idx:
				continue
			em.AddBond(newIndices[idx], newIndices[oidx], bond.GetBondType())
	res = em.GetMol()
	res.ClearComputedProps()
	Chem.GetSymmSSSR(res)
	res.UpdatePropertyCache(False)
	res._idxMap = newIndices
	return res


def _recursivelyModifyNs(mol, matches, indices=None):
	if indices is None:
		indices = []
	res = None
	while len(matches) and res is None:
		tIndices = indices[:]
		nextIdx = matches.pop(0)
		tIndices.append(nextIdx)
		nm = Chem.Mol(mol.ToBinary())
		nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
		nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
		cp = Chem.Mol(nm.ToBinary())
		try:
			Chem.SanitizeMol(cp)
		except ValueError:
			res, indices = _recursivelyModifyNs(nm, matches, indices=tIndices)
		else:
			indices = tIndices
			res = cp
	return res, indices


def AdjustAromaticNs(m, nitrogenPattern='[n&D2&H0;r5,r6]'):
	"""
		default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
		to fix: O=c1ccncc1
	"""
	try:
		Chem.GetSymmSSSR(m)
		m.UpdatePropertyCache(False)

		# break non-ring bonds linking rings:
		em = Chem.EditableMol(m)
		linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
		plsFix=set()
		for a,b in linkers:
			em.RemoveBond(a,b)
			plsFix.add(a)
			plsFix.add(b)
		nm = em.GetMol()
		for at in plsFix:
			at=nm.GetAtomWithIdx(at)
			if at.GetIsAromatic() and at.GetAtomicNum()==7:
				at.SetNumExplicitHs(1)
				at.SetNoImplicit(True)

		# build molecules from the fragments:
		fragLists = Chem.GetMolFrags(nm)
		frags = [FragIndicesToMol(nm,x) for x in fragLists]

		# loop through the fragments in turn and try to aromatize them:
		ok=True
		for i,frag in enumerate(frags):
			cp = Chem.Mol(frag.ToBinary())
			try:
				Chem.SanitizeMol(cp)
			except ValueError:
				matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
				lres,indices=_recursivelyModifyNs(frag,matches)
				if not lres:

					ok = False
					break
				else:
					revMap={}
					for k,v in frag._idxMap.iteritems():
						revMap[v]=k
					for idx in indices:
						oatom = m.GetAtomWithIdx(revMap[idx])
						oatom.SetNoImplicit(True)
						oatom.SetNumExplicitHs(1)
	except ValueError:
		ok = False

	if not ok:
		return None
	return m


def molPicker(m, printOut=False):
	unknMall = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
	if len(unknMall) == 1:

		numCs = len(unknMall[0].GetSubstructMatches(Cpatt))
		if numCs > 0:
			finalM = unknMall[0]
		else:
			finalM = None

	else:
		if printOut: print("  > %s molecules detected in the mol file!" % len(unknMall))

		iter2_Mols = []

		for unknMpotential in unknMall:
			unknMpotential.UpdatePropertyCache(strict=False)
			numCs = len(unknMpotential.GetSubstructMatches(Cpatt))

			if numCs > 0:
				iter2_Mols.append(unknMpotential)

		if len(iter2_Mols) == 1:
			if printOut: print("  > only 1 compound found containing carbon atoms")
			finalM = iter2_Mols[0]


		elif (len(iter2_Mols) < 3) & (len(iter2_Mols) > 1):  # if there is more than one molecule with carbon atoms in it, then resort to picking the one with the highest molecular weight
			if printOut: print("  > %s compounds found with carbon atoms, selecting primary compound by molecular weight" % len(iter2_Mols))
			iter3_Mols = []
			for unknMpotential in iter2_Mols:

				molwt = Descriptors.ExactMolWt(unknMpotential)
				iter3_Mols.append([molwt, unknMpotential])

			finalM = sorted(iter3_Mols, reverse=True)[0][1]

		else:
			if printOut: print("  <> %s organic compounds detected! Possible mixture, discarding..." % (len(iter2_Mols)))

			finalM = None

	return finalM


def superSanitizer(rawData, dataType="SMILES"):
	broken = False

	InChIKey = None
	SMILES = None

	if dataType == "SMILES":
		m = Chem.MolFromSmiles(rawData, sanitize=False)
	elif dataType == "MOL":
		m = Chem.MolFromMolBlock(rawData, sanitize=False)
	elif dataType == "SMARTS":
		m = Chem.MolFromSmarts(rawData)
	elif dataType == "rdkitMol":
		finalM = rawData

	if dataType != "rdkitMol":
		if m is not None:
			finalM = molPicker(m)
		else:
			finalM = None

	if finalM is not None:
		finalM.UpdatePropertyCache(strict=False)
		try:
			Chem.SanitizeMol(finalM)

			InChIKey = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(finalM)[0])
			SMILES = Chem.MolToSmiles(finalM)

		except ValueError:
			finalMcp = Chem.Mol(finalM.ToBinary())
			finalM_Nfixed = AdjustAromaticNs(finalMcp)

			if finalM_Nfixed is not None:
				finalM = finalM_Nfixed

				try:
					InChIKey = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(finalM)[0])
					SMILES = Chem.MolToSmiles(finalM)
					#print "\nN fixing worked!"
				except ValueError:
					pass

			if ((finalM_Nfixed is None) | (InChIKey is None)):
				#print "\nN fixing failed!"
				try:
					Chem.SanitizeMol(finalM,
						Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|
						Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
						Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=False)

					InChIKey = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(finalM)[0])
					SMILES = Chem.MolToSmiles(finalM)

				except ValueError as e:
					#print e
					#print "\nLast Ditch fixing!"
					try:
						Chem.SanitizeMol(finalM,
							Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|
							Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
							Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=False)

						finalM.UpdatePropertyCache(strict=False)
						Chem.SanitizeMol(finalM, Chem.SanitizeFlags.SANITIZE_SETAROMATICITY, catchErrors=True)

						InChIKey = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(finalM)[0])
						SMILES = Chem.MolToSmiles(finalM)

					except ValueError:
						broken = True

	return finalM, SMILES, InChIKey


def replacer(theM, theSites, internalTMS, nextIndex, theDepth, finalData, mode):
	#print "%s %s" % (theDepth, internalTMS-1)

	for i in range(nextIndex, len(theSites)):
		toReplace = theSites[i]

		SMARTSreplace = "[%s*]" % toReplace
		smartsM = Chem.MolFromSmarts(SMARTSreplace)

		replaceResult = Chem.ReplaceSubstructs(theM, smartsM, theTMS)[0]

		if theDepth < internalTMS-1:
			#print "loc 1"
			if mode == "Everything":
				aResult = Chem.Mol(replaceResult)  # you mant to make a copy of the molecule first

				for atom in aResult.GetAtoms():
					atom.SetIsotope(0)

				RR_wH = Chem.AddHs(aResult)
				RR_noH = Chem.RemoveHs(RR_wH)

				tmsInChIKey = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(RR_noH)[0])
				tmsSMILES = Chem.MolToSmiles(RR_noH)

				finalData.append((tmsInChIKey, tmsSMILES, theDepth+1))

			replacer(replaceResult, theSites, internalTMS, i+1, theDepth+1, finalData, mode)
			#print "finished an iteration"
		else:

			for atom in replaceResult.GetAtoms():
				atom.SetIsotope(0)

			RR_wH = Chem.AddHs(replaceResult)  # this seems redundant but I think it actually fixes some structural things
			RR_noH = Chem.RemoveHs(RR_wH)

			tmsInChIKey = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(RR_noH)[0])
			tmsSMILES = Chem.MolToSmiles(RR_noH)

			if mode == "Normal":
				finalData.append([tmsInChIKey, tmsSMILES, RR_wH])

			else: # partial!
				finalData.append((tmsInChIKey, tmsSMILES, internalTMS))


if __name__ == '__main__':

	with open(sys.argv[1]) as f:

		for line in f:
			rawSMILES = line.rstrip()

			subStructInChIKeys = []

			finalM, finalSMILES, finalInChIKey = superSanitizer(rawSMILES, "SMILES")

			print("Raw Input SMILES: %s Computed InChIKey: %s" % (rawSMILES, finalInChIKey))

			finalM_wH = Chem.AddHs(finalM)

			OHmatches = finalM_wH.GetSubstructMatches(theOH)

			#if derivedPats[theSelection] == '5523':
			#	print Chem.MolToSmiles(finalM_wH)

			if len(OHmatches) > 0:
				Hlabels_OH = []
				for aMatch in OHmatches:
					idxH = aMatch[1]  # this should be the index of the hydrogren atom that you'll be replacing
					Hatom = finalM_wH.GetAtomWithIdx(idxH)
					Hatom.SetIsotope(idxH)

					Hlabels_OH.append(idxH)

				finalData = []
				replacer(finalM_wH, Hlabels_OH, len(OHmatches), 0, 0, finalData, "Normal")

				#print finalData[0]
				subStructInChIKeys.append((finalData[0][0], finalData[0][1], len(OHmatches)))
				finalM_wH = finalData[0][2]

			SHmatches = finalM_wH.GetSubstructMatches(theSH)
			NHmatches = finalM_wH.GetSubstructMatches(theNH)
			Hlabels_SandN = []

			for aMatch in SHmatches:
				idxH = aMatch[1]  # this should be the index of the hydrogren atom that you'll be replacing
				Hatom = finalM_wH.GetAtomWithIdx(idxH)
				Hatom.SetIsotope(idxH)

				Hlabels_SandN.append(idxH)

			for aMatch in NHmatches:
				idxH = aMatch[1]  # this should be the index of the hydrogren atom that you'll be replacing
				Hatom = finalM_wH.GetAtomWithIdx(idxH)
				Hatom.SetIsotope(idxH)

				Hlabels_SandN.append(idxH)

			if len(Hlabels_SandN) > 0:
				finalData = []
				replacer(finalM_wH, Hlabels_SandN, len(Hlabels_SandN), 0, 0, finalData, "Everything")
				subStructInChIKeys.extend(finalData)

			totalSum = len(OHmatches) + len(SHmatches) + len(NHmatches)

			print(" >> There are %s OH sites, %s SH sites, and %s NH sites for TMS derivatization, %s potential derivatizations" % (len(OHmatches), len(SHmatches), len(NHmatches), len(subStructInChIKeys)))

			uniqStructures = dict()
			for derivInChIKey, derivSMILES, derivTMScount in subStructInChIKeys:
				shortKey = derivInChIKey.split("-")[0]
				if shortKey not in uniqStructures:
					uniqStructures[shortKey] = [derivInChIKey, derivSMILES, derivTMScount]

					print("   %s TMS groups, InChIKey: %s\tSMILES: %s" % (derivTMScount, derivInChIKey, derivSMILES))


			print(" >> %s unique TMS derivatives found\n" % (len(uniqStructures)))
