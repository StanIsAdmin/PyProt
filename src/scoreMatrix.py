from sequence import Sequence
from aminoAcid import AminoAcid

class ScoreMatrix:
	"""
	Represents a scoring matrix, used to determine the score between two Amino Acids
	"""
	
	def __init__(self, path="", description="", ignore=""):
		"""
		Creates a ScoreMatrix object.
		If 'path' is provided, loads the Score values from an iij file.
		Otherwise, creates a ScoreMatrix for all possible AminoAcids with values 0.
		"""
		self._description = description
		self._ignore = Sequence(ignore)
		self._matrix = []
		self._aaOrder = {}
		self._aaSequence = Sequence()
		
		#If path is provided, load directly from iij file
		if path != "":
			with open(path, 'r') as file:
				foundAAOrder = False #Have we found the line with the amino acid values and order yet?
				for line in file:
					if line[0] != "#": #Comments
						
						if not foundAAOrder: #Read aa values and order
							for aa in line.split():
								self._aaSequence.extend(aa)
							self._aaOrder = {aa: index for aa, index in zip(self._aaSequence, range(len(self._aaSequence)))}
							foundAAOrder = True
						else: #Read matrix values
							self._matrix.append([int(v) for v in line.split()])
		
		#Otherwise initialize matrix with 0
		else:
			lineSize = 1
			for aa in AminoAcid.getAllNames():
				if AminoAcid(aa) not in self._ignore:
					self._aaSequence.extend(aa)
					self._aaOrder[self._aaSequence[-1]] = lineSize-1
					self._matrix.append([0 for i in range(lineSize)])
					lineSize += 1
				
	#Representation
	def __repr__(self):
		sepSize = 4
		result = [" - - - " + self._description + " - - - "]
		for values, aa in zip(self._matrix, self._aaSequence):
			tempstr = '{a!s:<{w}}'.format(a=aa, w=sepSize)
			for value in values:
				tempstr += '{v:<{w}}'.format(v=value, w=sepSize)
			result.append(tempstr)
		tempstr = " "*sepSize
		for aa in self._aaSequence :
			tempstr += '{a!s:<{w}}'.format(a=aa, w=sepSize)
		result.append("")
		result.append(tempstr)
		return "\n".join(result)
	
	#Scoring
	def setScore(self, aa1, aa2, score):
		id1 = self._aaOrder[aa1]
		id2 = self._aaOrder[aa2]
		if id1 > id2:
			self._matrix[id1][id2] = score
		else:
			self._matrix[id2][id1] = score
		
	def getScore(self, aa1, aa2):
		id1 = self._aaOrder[aa1]
		id2 = self._aaOrder[aa2]
		if id1 > id2:
			return self._matrix[id1][id2]
		else:
			return self._matrix[id2][id1]

