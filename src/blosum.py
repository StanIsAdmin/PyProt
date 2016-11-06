from aminoAcid import AminoAcid
from sequence import Sequence, loadFasta
from scoreMatrix import ScoreMatrix
from alignMatrix import AlignMatrix
from math import log, ceil
import cProfile

def belongs(outSequence, group, minMatches):
	"""
	Returns True if there's at least 'minMatches' matches between outSequence and any sequence 
	from 'group'; False otherwise.
	Size of all sequences is assumed to be equal.
	"""
	for inSequence in group:
		matches = 0
		for i in range(len(outSequence)):
			if outSequence[i] == inSequence[i]:
				matches += 1
			if matches >= minMatches: #As soon as a match is found
				return True
	return False

	
def makeGroupsFromFasta(path, requiredIdentityPercent):
	"""
	Loads Sequences from file at 'path' and separate them in groups with an identity of 
	at least 'requiredIdentityPercent'.
	Returns a list representing the groups as lists of Sequence objects.
	"""
	sequences = [seq for seq in loadFasta(path)]
	groups = [[sequences[0]]] #First sequence is assigned to first group
	
	seqSize = len(sequences[0]) #Size of the sequences
	
	#Number of matches required to achieve requiredIdentityPercent
	minMatches = ceil((requiredIdentityPercent/100)*seqSize)
	
	#For each outSequence not yet in a group
	for outSequence in sequences[1:]:
		groupFound = False #has a group been found ?
		groupIndex = 0 #index of the group we're looking in
		
		#Look for a group where outSequence belongs
		while (not groupFound) and groupIndex < len(groups):
			if belongs(outSequence, groups[groupIndex], minMatches):
				groups[groupIndex].append(outSequence)
				groupFound = True
			else:
				groupIndex += 1 #Move on to next group
		
		#If no group works, create a new one
		if not groupFound:
			groups.append([outSequence])
	
	return groups

	
def valueDictsFromGroups(groups):
	"""
	Transforms each group from 'groups' into a list of dictionaries, one for each Sequence column.
	The dictionaries map each AminoAcid found in that column to their count.
	Returns the list of dictionaries and a list of the size (in Sequences) of their groups.
	"""
	seqSize = len(groups[0][0])
	groupCount = len(groups)
	
	groupValues = []
	groupSizes = []
	
	for group in groups: #For each group
		groupAAs = []
		groupSize = len(group) #Size of group (n°sequences)
		
		for col in range(seqSize): #For each column
			groupCol = {}
			
			for seq in group: #For each sequence in group
				try:
					groupCol[seq[col]] += 1
				except:
					groupCol[seq[col]] = 1
			
			groupAAs.append(groupCol)
		groupValues.append(groupAAs)
		groupSizes.append(groupSize)
	return groupValues, groupSizes

	
def getFrequencies(groupValues, groupSizes):
	seqSize = len(groupValues[0])
	groupCount = len(groupSizes)
	freqPairs = {} #frequencies of amino acid pairs (fAB)
	
	#Frequencies of single amino acids (fA)
	freqSingle = {AminoAcid(aa):0 for aa in AminoAcid.getAllNames()} 
	
	freqSum = 0 #Sum of frequencies  sum(fAB)
	
	for col in range(seqSize):
		for groupAIndex in range(groupCount-1):
			groupA = groupValues[groupAIndex]
			groupASize = groupSizes[groupAIndex]
			
			for groupBIndex in range(groupAIndex+1, groupCount):
				groupB = groupValues[groupBIndex]
				groupBSize = groupSizes[groupBIndex]
				
				for aaA, aaACount in groupA[col].items():
					aaAPercent = aaACount / groupASize
					
					for aaB, aaBCount in groupB[col].items():
						aaBPercent = aaBCount / groupBSize
						aaScore = aaAPercent * aaBPercent
						
						freqSum += aaScore	#Sum of frequencies			
						freqSingle[aaA] += aaScore/2
						freqSingle[aaB] += aaScore/2
						
						pairIndex = (aaA, aaB) if aaA > aaB else (aaB, aaA)
						try:
							freqPairs[pairIndex] += aaScore
						except:
							freqPairs[pairIndex] = aaScore
	return freqPairs, freqSingle, freqSum
	
def sumFrequenciesToProb(freqPairsList, freqSingleList, freqSumList):
	fSum = sum(freqSumList)
	
	probPairs, probSingle = {}, {}
	for freqPairs, freqSingle in zip(freqPairsList, freqSingleList):
		for key,value in freqPairs.items():
			try:
				probPairs[key] += value / fSum
			except:
				probPairs[key] = value / fSum
		for key,value in freqSingle.items():
			try:
				probSingle[key] += value / fSum
			except:
				probSingle[key] = value / fSum
	return probPairs, probSingle
		
	
def blosumFromFrequencies(scoreMatrix, probPairs, probSingle):
	for key, qAB in probPairs.items():
		
		aaA, aaB = key #Both amino acids
		pA, pB = probSingle[aaA],probSingle[aaB] #Their single probabilities
		
		eAB = pA * pA if aaA == aaB else 2 * pA * pB
		sAB = int(round(2 * log(qAB / eAB, 2)))
		scoreMatrix.setScore(aaA, aaB, sAB)
	
	return scoreMatrix
			
		
			
def blosumFromFasta(requiredIdentityPercent, *filepaths):
	freqPairsList, freqSingleList, freqSumList = [], [], []
	#MULTITHREAD
	for path in filepaths:
		groups = makeGroupsFromFasta(path, requiredIdentityPercent)
		groupValues,groupSizes = valueDictsFromGroups(groups)
		freqPairs, freqSingle, freqSum = getFrequencies(groupValues, groupSizes)
		
		freqPairsList.append(freqPairs) #Not atomic, only addition later
		freqSingleList.append(freqSingle)
		freqSumList.append(freqSum)
	#END MULTITHREAD
	probPairs, probSingle = sumFrequenciesToProb(freqPairsList, freqSingleList, freqSumList)
	scoreMatrix = ScoreMatrix("", "BLOSUM{}".format(requiredIdentityPercent), "BZJUO")
	blosum = blosumFromFrequencies(scoreMatrix, probPairs, probSingle)
	print(blosum)

	
path1 = r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\SH3-A.fasta"
path2 = r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\SH3-B.fasta"
path3 = r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\SH3-C.fasta"
path4 = r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\SH3-D.fasta"
blosumFromFasta(70, path1)