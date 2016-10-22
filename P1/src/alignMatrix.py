from scoreMatrix import ScoreMatrix
from sequence import Sequence
from aminoAcid import AminoAcid

class AlignMatrix:
	"""
	Represents an alignment matrix, used to determine the alignment score between two sequences
	"""
	
	def __init__(self, scoreMatrix, seqA, seqB):
		if not (isinstance(seqA, Sequence) and isinstance(seqB, Sequence)):
			raise TypeError("seqA and seqB must be Sequence objects")
		if len(seqA)==0 or len(seqB)==0:
			raise ValueError("seqA and seqB cannot be empty")
			
		if not isinstance(scoreMatrix, ScoreMatrix):
			raise TypeError("scoreMatrix must be a Score object")
			
		self._rowSeq = seqA
		self._colSeq = seqB
		self._scoreMatrix = scoreMatrix
		
		self._gapPenalty = -4
		self._alignMatrix = [[self._gapPenalty*i for i in range(len(self._colSeq)+1)]]
		self._alignMatrix.extend([[(self._gapPenalty*j if i==0 else 0) for i in range(len(self._colSeq)+1)] for j in range(1, len(self._rowSeq))])
		
	def fill(self):
		pass
		
		
	def __repr__(self):
		sepSize = 5
		print(" "*2*sepSize, end="")
		for aa in self._colSeq :
			print('{v!s:<{w}}'.format(v=aa, w=sepSize), end="")
		print("\n", end=" "*sepSize)
		for value in self._alignMatrix[0]:
			print('{v:<{w}}'.format(v=value, w=sepSize), end="")
		print()
		for line, aa in zip(self._alignMatrix[1:], self._rowSeq):
			print(aa, end=" "*(sepSize-1))
			for value in line:
				print('{v:<{w}}'.format(v=value, w=sepSize), end="")
			print()
		return ""
		

s = ScoreMatrix(r"..\Resources\blosum\blosum30.iij")
print(s)
a = AlignMatrix(s, Sequence("shortGGGTTF"), Sequence("shortMGGETFA"))
print(a)
		