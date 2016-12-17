from math import log
from aminoacid import AminoAcid


class GOR3:
	"""
	Implements the GOR III secondary structure prediction algorithm.
	Objects of this class must be trained with known examples of sequences and their structure
	before being able to predict the structures of new sequences.
	"""

	def __init__(self):
		self.structures = "HETC"
		self.neighbourOffset = 8
		
		self.trainings = 0 #Number of trainings (one per AA)
		self.strucCount = {s:0 for s in self.structures}
		self.pairCount = {(s, a):0 for s in self.structures for a in AminoAcid.getAllNames()}
		self.tripletCount = {}
		for s in self.structures:
			for a in AminoAcid.getAllNames():
				for na in AminoAcid.getAllNames(): #Neighbour AA
					self.tripletCount[(s, a, na)] = 0
		
		
	def train(self, sequence, structure):
		"""
		Trains the system with a known example of a sequence and its structure.
		"""
		self.trainings += len(sequence)
		
		for index in range(len(sequence)):
			curAminoacid = sequence[index]
			curStructure = structure[index]
			self.strucCount[curStructure] += 1
			self.pairCount[(curStructure, curAminoacid)] += 1
			
			for neiAminoacid in self.neighbourValues(sequence, index):
				self.tripletCount[(curStructure, curAminoacid, neiAminoacid)] += 1
	
	
	def predict(self, sequence):
		"""
		Returns the predicted structure of 'sequence', based on received training.
		"""
		structure = [] #Result: predicted structure
		
		#Predict structures for each aminoacid in sequence
		for index in range(len(sequence)):
			curAminoacid = sequence[index]
			
			#First possible structure
			predStructure = self.structures[0]
			predScore = self.__getScore(sequence, index, predStructure)
			
			#Other structures
			for curStructure in self.structures[1:]:
				curScore = self.__getScore(sequence, index, curStructure)
				
				#Remember structure that gives best score
				if curScore > predScore:
					predStructure = curStructure
					predScore = curScore
					
			structure.append(predStructure)
		return "".join(structure)
		
	
	def __getScore(self, sequence, index, struct):
		"""
		Returns I(deltaS, R) as defined by the GOR III algorithm.
		"""
		aminoacid = sequence[index]
		scoreTerms = []
		
		score = self.pairCount[(struct, aminoacid)]
		score /= sum(self.pairCount[(otherStruct, aminoacid)] for otherStruct in self.getStructures(struct))
		scoreTerms.append(score)
		
		score = (self.trainings - self.strucCount[struct]) / self.strucCount[struct]
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