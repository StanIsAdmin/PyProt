from sequence import Sequence
from aminoAcid import AminoAcid

class ScoreMatrix:
	"""
	Represents a scoring matrix, used to determine the score between two Amino Acids
	"""
	
	def __init__(self, path=""):
		self._matrix = []
		self._aaOrder = {}
		self._aaValues = []
		
		#If path is provided, load directly from file
		if path != "":
			with open(path, 'r') as file:
				foundAAOrder = False
				for line in file:
					if line[0] != "#":
						if not foundAAOrder:
							self._aaValues = line.split()
							self._aaOrder = {value: index for value, index in zip(self._aaValues, range(len(self._aaValues)))}
							foundAAOrder = True
						else:
							self._matrix.append([int(v) for v in line.split()])
		#Otherwise initialize with 0
		else:
			lineSize = 1
			for aa in AminoAcid.getAllNames():
				self._matrix.append([0 for i in range(lineSize)])
				lineSize += 1
				
	def __repr__(self):
		print(self._aaValues)
		for line in self._matrix:
			print(line)
		return 
		
	def __repr__(self):
		sepSize = 4
		
		for values, aa in zip(self._matrix, self._aaValues):
			print('{v!s:<{w}}'.format(v=aa, w=sepSize), end="")
			for value in values:
				print('{v:<{w}}'.format(v=value, w=sepSize), end="")
			print()
		print(" "*sepSize, end="")
		for aa in self._aaValues :
			print('{v!s:<{w}}'.format(v=aa, w=sepSize), end="")
		
		return ""
			
	def setScore(self, score, aa1, aa2):
		id1 = self._aaOrder[str(aa1)]
		id2 = self._aaOrder[str(aa2)]
		if id1 > id2:
			self._matrix[id1][id2] = score
		else:
			self._matrix[id2][id1] = score
		
	def getScore(self, aa1, aa2):
		id1 = self._aaOrder[str(aa1)]
		id2 = self._aaOrder[str(aa2)]
		if id1 > id2:
			return self._matrix[id1][id2]
		else:
			return self._matrix[id2][id1]

				
		
			
"""
ScoreMatrix(r"Resources\blosum\blosum30.iij")
"""