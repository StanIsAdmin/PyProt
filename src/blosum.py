from sequence import Sequence, loadFasta
from scoreMatrix import ScoreMatrix
from alignMatrix import AlignMatrix
from math import ceil


def makeGroupsFromFasta(path, requiredIdentityPercent):
	sequences = [seq for seq in loadFasta(path)]
	groups = [[sequences[0]]]
	
	#For each outSequence not yet in a group
	for outSequence in sequences[1:]:
		groupFound = False #has a group been found ?
		groupIndex = 0 #index of the group we're looking in
		
		#Look for a group where outSequence belongs
		while (not groupFound) and groupIndex < len(groups):
			if belongs(outSequence, groups[groupIndex], requiredIdentityPercent):
				groups[groupIndex].append(outSequence)
				groupFound = True
			else:
				groupIndex += 1 #Move on to next group
		
		#If no group works, create a new one
		if not groupFound:
			groups.append([outSequence])
	
	return groups

def belongs(outSequence, group, requiredIdentityPercent):
	requiredIdentity = requiredIdentityPercent/100

	for inSequence in group:
		minSize = min(len(outSequence), len(inSequence))
		minMatches = ceil(requiredIdentity*minSize)
		matches = 0
		for i in range(minSize):
			if outSequence[i] == inSequence[i]:
				matches += 1
			if matches >= minMatches:
				return True
	return False
			

path = r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\test.fasta"
groups = makeGroupsFromFasta(path, 50)
for group in groups:
	print(group)