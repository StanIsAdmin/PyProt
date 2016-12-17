from math import log
from aminoacid import AminoAcid


class GOR3:
	def __init__(self):
		self.structures = "HETC"
		self.neighbourOffset = 8
		
		self.trainings = 0
		self.structureCount = {s:0 for s in self.structures}
		self.pairCount = {(s, a):0 for s in self.structures for a in AminoAcid.getAllNames()}
		self.tripletCount = {}
		for s in self.structures:
			for a in AminoAcid.getAllNames():
				for na in AminoAcid.getAllNames(): #Neighbour AA
					self.tripletCount[(s, a, na)] = 0
		
		
	def train(self, sequence, structure):
		self.trainings += len(sequence)
		for index in range(len(sequence)):
			curAminoacid = sequence[index]
			curStructure = structure[index]
			self.structureCount[curStructure] += 1
			self.pairCount[(curStructure, curAminoacid)] += 1
			
			for neiAminoacid in self.neighbourValues(sequence, index):
				self.tripletCount[(curStructure, curAminoacid, neiAminoacid)] += 1
	
	
	def predict(self, sequence):
		structure = [] #Result: predicted structure
		
		for index in range(len(sequence)):
			curAminoacid = sequence[index]
			
			predStructure = ""
			predScore = -1
			for curStructure in self.structures:
				curScore = self.__getScore(sequence, index, curStructure)
				
				if curScore > predScore:
					predStructure = curStructure
					predScore = curScore
					
			structure.append(predStructure)
		return "".join(structure)
		
	
	def __getScore(self, sequence, index, struct):
		aminoacid = sequence[index]
		scoreTerms = []
		
		score = self.pairCount[(struct, aminoacid)]
		score /= sum(self.pairCount[(otherStruct, aminoacid)] for otherStruct in self.getStructures(struct))
		scoreTerms.append(score)
		
		score = (self.trainings - self.structureCount[struct]) / self.structureCount[struct]
		scoreTerms.append(score)
		
		for neiAminoacid in self.neighbourValues(sequence, index):
			score = self.tripletCount[(struct, aminoacid, neiAminoacid)]
			score /= sum(self.tripletCount[(otherStruct, aminoacid, neiAminoacid)] for otherStruct in self.getStructures(struct))
			scoreTerms.append(score)
			
			score = sum(self.pairCount[(otherStruct, aminoacid)] for otherStruct in self.getStructures(struct))
			score /= self.pairCount[(struct, aminoacid)]
			scoreTerms.append(score)
			
		return sum(log(s) for s in scoreTerms)
			
			
		
	def neighbourOffsets(self):
		for offset in range(-self.neighbourOffset, self.neighbourOffset+1):
			if offset!=0:
				yield offset
	
	
	def neighbourValues(self, sequence, index):
		for offset in self.neighbourOffsets():
			neiIndex = index + offset
			if neiIndex >= 0 and neiIndex < len(sequence):
				yield sequence[neiIndex]
	
	
	def getStructures(self, exclude=None):
		for s in self.structures:
			if s != exclude:
				yield s