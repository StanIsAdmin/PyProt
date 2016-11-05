from sequence import Sequence, loadFasta
from scoreMatrix import ScoreMatrix
from alignMatrix import AlignMatrix
from math import ceil


"""
def makeGroupsFromFasta(path, requiredIdentity):
	sequences = [seq for seq in loadFasta(path)]
	for seq in sequences:
		print(seq.getDescription())
		print(seq)
	
	
	groups = [[sequences[0]]]
	return groups
	for seqA in sequences[1:]:
		
		groupFound = False
		groupIndex = 0
		group = groups[groupIndex]
		
		while (not groupFound) and len(groups) > groupIndex:
			belongs = True
			seqIndex = 0
			seqB = group[seqIndex]
			while belongs and len(group) > seqIndex:
				alignment = align.firstGlobalAlign(seqA,seqB,1)
				if alignment.gaps == 0 and alignment.identityPercent >= requiredIdentity:
					seqIndex += 1
					seqB = group[seqIndex]
				else:
					belongs = False
			if belongs:
				group.append(seqA)
				groupFound = True
				groupIndex = 0
			else:
				groupIndex += 1
			
		if not groupFound:
			groups.append([seqA])
	
	return groups
"""
def makeGroupsFromFasta(path, requiredIdentity):
	sequences = [seq for seq in loadFasta(path)]
	groups = [[sequences[0]]]
	align = AlignMatrix(ScoreMatrix())
	
	#For each outSequence not yet in a group
	for outSequence in sequences[1:]:
		groupFound = False #has a group been found ?
		groupIndex = 0 #index of the group we're looking in
		
		#Look for a group where outSequence belongs
		while (not groupFound) and groupIndex < len(groups):
			group = groups[groupIndex] #group we're looking in
			belongs = False #does the sequence belong in this group ?
			seqIndex = 0 #index of the sequence we're comparing it to
			
			while (not belongs) and seqIndex < len(group):
				inSequence = group[seqIndex] #sequence from the group
				
				#Check if any sequences "match" the outSequence
				alignment = align.firstGlobalAlign(outSequence,inSequence,-1000)
				if alignment.gaps == 0 and alignment.identityPercent >= requiredIdentity:
					belongs = True
				else:
					print(alignment)
					seqIndex += 1 #Move on to next sequence
			if belongs:
				group.append(outSequence)
				groupFound = True
			else:
				groupIndex += 1 #Move on to next group
		
		#If no group works, create a new one
		if not groupFound:
			groups.append([outSequence])
	
	return groups
	


path = r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\test.fasta"
groups = makeGroupsFromFasta(path, 50)
for group in groups:
	print(group)