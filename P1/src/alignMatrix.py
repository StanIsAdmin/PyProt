from scoreMatrix import BlosumMatrix
from sequence import Sequence
from aminoAcid import AminoAcid

class AlignMatrix:
	"""
	Represents an alignment matrix, used to determine the alignment score between two sequences
	"""
	
	def __init__(self, blosumMatrix, seqA, seqB):
		if not (isinstance(seqA, Sequence) and isinstance(seqB, Sequence)):
			raise TypeError("seqA and seqB must be Sequence objects")
		if len(seqA)==0 or len(seqB)==0:
			raise ValueError("seqA and seqB cannot be empty")
			
		if not isinstance(blosumMatrix, BlosumMatrix):
			raise TypeError("scoreMatrix must be a Score object")
			
		self._rowSeq = seqA
		self._colSeq = seqB
		self._blosumMatrix = blosumMatrix
		
		self._gapPenalty = -4
		self._alignMatrix = [[self._gapPenalty*i for i in range(len(self._colSeq)+1)]]
		self._alignMatrix.extend([[(self._gapPenalty*j if i==0 else 0) for i in range(len(self._colSeq)+1)] for j in range(1, len(self._rowSeq)+1)])
		
	def fill(self):
		for row in range(1, len(self._rowSeq)+1):
			for col in range(1, len(self._colSeq)+1):
				v1 = self._alignMatrix[row-1][col] + self._gapPenalty
				v2 = self._alignMatrix[row][col-1] + self._gapPenalty
				v3 = self._alignMatrix[row-1][col-1] + \
					self._blosumMatrix.getScore(self._rowSeq[row-1], self._colSeq[col-1])
				self._alignMatrix[row][col] = max(v1, v2, v3)
		
		
	def __repr__(self):
		sepSize = 5
		#Amino Acid column sequence
		print(" "*2*sepSize, end="")
		for aa in self._colSeq :
			print('{v!s:<{w}}'.format(v=aa, w=sepSize), end="")
		print("\n", end=" "*sepSize)
		#First line of values (no line sequence)
		for value in self._alignMatrix[0]:
			print('{v:<{w}}'.format(v=value, w=sepSize), end="")
		print()
		#Amino Acid line sequence and values
		for line, aa in zip(self._alignMatrix[1:], self._rowSeq):
			print('{v!s:<{w}}'.format(v=aa, w=sepSize), end="")
			for value in line:
				print('{v:<{w}}'.format(v=value, w=sepSize), end="")
			print()
		return ""
		

b = BlosumMatrix(r"..\Resources\blosum\blosum62.iij")
print(b)
a = AlignMatrix(b, Sequence("shortGGVTTF"), Sequence("shortMGGETFA"))
print(a)
a.fill()
print(a)
		