from sequence import Sequence
from aminoAcid import AminoAcid

class ScoreMatrix:
	"""
	Represents a scoring matrix, used to determine the score between two Amino Acids
	"""
	
	def __init__(self, path=""):
		self._matrix = []
		
		#If path is provided, load directly from file
		if path != "":
			with open(path, 'r') as file:
				for line in file:
					if line[0] != "#" and line.strip()[0] != "A":
						self._matrix.append([int(v) for v in line.split()])
		#Otherwise initialize with 0
		else:
			lineSize = 1
			for aa in AminoAcid.getAllNames():
				self._matrix.append([0 for i in range(lineSize)])
				lineSize += 1
				
	def __repr__(self):
		for line in self._matrix:
			print(line)
		return ""
			
	def setScore(self, score, aa1, aa2):
		id1 = aa1.getId() if isinstance(aa1, AminoAcid) else AminoAcid.getIdByName(aa1)
		id2 = aa2.getId() if isinstance(aa2, AminoAcid) else AminoAcid.getIdByName(aa2)
		if id1 > id2:
			self._matrix[id1][id2] = score
		else:
			self._matrix[id2][id1] = score
		
	def getScore(self, aa1, aa2):
		id1 = aa1.getId() if isinstance(aa1, AminoAcid) else AminoAcid.getIdByName(aa1)
		id2 = aa2.getId() if isinstance(aa2, AminoAcid) else AminoAcid.getIdByName(aa2)
		if id1 > id2:
			return self._matrix[id1][id2]
		else:
			return self._matrix[id2][id1]

				
		
			
"""
ScoreMatrix(r"Resources\blosum\blosum30.iij")
"""