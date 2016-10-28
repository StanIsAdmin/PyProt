from scoreMatrix import ScoreMatrix
from sequence import Sequence
from aminoAcid import AminoAcid
from copy import deepcopy

class AlignMatrix:
	"""
	Represents an alignment matrix, used to determine the alignment score between two sequences
	"""
	
	def __init__(self, scoreMatrix, seqA, seqB, initGapPenalty, extGapPenalty):
		if not (isinstance(seqA, Sequence) and isinstance(seqB, Sequence)):
			raise TypeError("seqA and seqB must be Sequence objects")
		if len(seqA)==0 or len(seqB)==0:
			raise ValueError("seqA and seqB cannot be empty")
			
		if not isinstance(scoreMatrix, ScoreMatrix):
			raise TypeError("scoreMatrix must be a Score object")
			
		self._rowSeq = seqA
		self._colSeq = seqB
		self._scoreMatrix = scoreMatrix
		
		self._maxAlignRow = 0
		self._maxAlignCol = 0
		
		self._initGapPenalty = initGapPenalty
		self._extGapPenalty = extGapPenalty
		
		self._alignMatrix = [[0]]
		self._alignMatrix[0].extend([self._initGapPenalty+self._extGapPenalty*(i-1) for i in range(1, len(self._colSeq)+1)])
		self._alignMatrix.extend([[(self._initGapPenalty+self._extGapPenalty*(j-1) if i==0 else 0) for i in range(len(self._colSeq)+1)] for j in range(1, len(self._rowSeq)+1)])
		
		self._rowGapMatrix = deepcopy(self._alignMatrix)
		self._colGapMatrix = deepcopy(self._alignMatrix)
	
	
	def __fill(self, local):
		self._maxAlignRow = 0
		self._maxAlignCol = 0
		for row in range(1, len(self._rowSeq)+1):
			for col in range(1, len(self._colSeq)+1):				
				#Row gap matrix
				initialGap = self._alignMatrix[row-1][col] + self._initGapPenalty
				extendedGap = self._rowGapMatrix[row-1][col] + self._extGapPenalty
				self._rowGapMatrix[row][col] = max(initialGap, extendedGap)
				
				#Col gap matrix
				initialGap = self._alignMatrix[row][col-1] + self._initGapPenalty
				extendedGap = self._colGapMatrix[row][col-1] + self._extGapPenalty
				self._colGapMatrix[row][col] = max(initialGap, extendedGap)
				
				#Align matrix
				insertScore = self._rowGapMatrix[row][col]
				deleteScore = self._colGapMatrix[row][col]
				matchScore = self._alignMatrix[row-1][col-1] + \
					self._scoreMatrix.getScore(self._rowSeq[row-1], self._colSeq[col-1])
					
				if local:
					self._alignMatrix[row][col] = max(insertScore, deleteScore, matchScore, 0)
				else:
					self._alignMatrix[row][col] = max(insertScore, deleteScore, matchScore)
				
				if self._alignMatrix[row][col] > self._alignMatrix[self._maxAlignRow][self._maxAlignCol]:
					self._maxAlignRow = row
					self._maxAlignCol = col
		print(self)
	
	
	def __align(self, i, j, seqA, seqB):
		
		if i>0 or j>0:
			
			if i > 0 and j > 0 and self._alignMatrix[i][j] == \
			self._alignMatrix[i-1][j-1] + self._scoreMatrix.getScore(self._rowSeq[i-1], self._colSeq[j-1]):
				seqA.insert(self._rowSeq[i-1], 0)
				seqB.insert(self._colSeq[j-1], 0)
				yield from self.__align(i-1, j-1, seqA, seqB)
				seqA.delete(0)
				seqB.delete(0)
			
			if i > 0 and self._alignMatrix[i][j] == self._alignMatrix[i-1][j] + self._initGapPenalty:
				seqA.insert(self._rowSeq[i-1], 0)
				seqB.insert("|", 0)
				yield from self.__align(i-1, j, seqA, seqB)
				seqA.delete(0)
				seqB.delete(0)
				
			if j > 0 and self._alignMatrix[i][j] == self._alignMatrix[i][j-1] + self._initGapPenalty:
				seqA.insert("|", 0)
				seqB.insert(self._colSeq[j-1], 0)
				yield from self.__align(i, j-1, seqA, seqB)
				seqA.delete(0)
				seqB.delete(0)
		else:
			yield Sequence(seqA), Sequence(seqB)
			
	def globalAlign(self):
		self.__fill(False) #Fill in global mode
		yield from self.__align(len(self._rowSeq), len(self._colSeq), Sequence(), Sequence())
	
	def localAlign(self):
		self.__fill(True) #Fill in local mode
		yield from self.__align(self._maxAlignRow, self._maxAlignCol, Sequence(), Sequence())
		
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
		

s = ScoreMatrix(r"..\Resources\blosum\blosum62.iij")
print(s)
a = AlignMatrix(s, Sequence("shortGGVTTF"), Sequence("shortMGGETFA"), -4, -1)
print(a)

for seqA, seqB in a.localAlign():
	print(seqA, seqB, sep="\n")
		