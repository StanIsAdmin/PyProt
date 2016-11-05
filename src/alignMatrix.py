from scoreMatrix import ScoreMatrix
from sequence import Sequence
from aminoAcid import AminoAcid
from copy import deepcopy

class AlignedSequences:
	def __init__(self, seqA, seqB, score):
		self.seqA = seqA
		self.seqB = seqB
		self.score = score

class AlignMatrix:
	"""
	Represents an alignment matrix, used to determine the alignment score between two sequences
	"""
	
	def __init__(self, scoreMatrix, seqA, seqB):
		"""
		Creates an AlignMatrix object that uses scoreMatrix as a scoring system between amino acids from sequences seqA and seqB.
		"""
		if not (isinstance(seqA, Sequence) and isinstance(seqB, Sequence)):
			raise TypeError("seqA and seqB must be Sequence objects")
		if len(seqA)==0 or len(seqB)==0:
			raise ValueError("seqA and seqB cannot be empty")
			
		if not isinstance(scoreMatrix, ScoreMatrix):
			raise TypeError("scoreMatrix must be a ScoreMatrix object")
			
		self._colSeq = Sequence(seqA) #copies of sequences to align
		self._rowSeq = Sequence(seqB)
		self._scoreMatrix = scoreMatrix #scoring matrix to use
		
		#Used for alignment scores between sequences
		self._alignMatrix = [[0 for i in range(len(self._colSeq)+1)] for j in range(len(self._rowSeq)+1)]
		
		#same kinds of matrix as _alignMatrix, used for affine gap penalty
		self._rowGapMatrix = deepcopy(self._alignMatrix)
		self._colGapMatrix = deepcopy(self._alignMatrix)
		
		#Used for backtracking purposes
		self._originMatrix = [["" for i in range(len(self._colSeq)+1)] for j in range(len(self._rowSeq)+1)]
		
		self._alignMode = "" #Align mode (global, semiglobal, local)
		
		self._maxScoreRow = 0 #Indexes of maximum value from _alignMatrix
		self._maxScoreCol = 0
		
		self._maxAlignScore = 0
		self._currentAlignPath = []
		self._bestAlignPath = [] #sequence of origins for alignment with best score
		
		self._iniGapPenalty = 0 #Initial gap penalty
		self._extGapPenalty = 0 #Extended gap penalty
		
		#Aligned sequences (result of alignment)
		self._alignedSeqA = None
		self._alignedSeqB = None
		
		
	def globalAlign(self, iniGapPenalty, extGapPenalty=None, semiGlobal=False):
		if semiGlobal:
			self._alignMode = "semiglobal"
		else:
			self._alignMode = "global"
			
		self.__initialize(iniGapPenalty, extGapPenalty) #Initialize all data structures
		print(self)
		#Align in global mode, from the bottom right
		yield from self.__align(len(self._rowSeq), len(self._colSeq))
	
	
	def localAlign(self, iniGapPenalty, extGapPenalty=None, subOptimal=False):
		self._alignMode = "local"
		self.__initialize(iniGapPenalty, extGapPenalty) #Initialize all data structures
		print(self)
		#Align in local mode, from the maximum value (optimal)
		yield from self.__align(self._maxScoreRow, self._maxScoreCol)
		
		#If suboptimal alignments are requested
		if subOptimal:
			self.__clearBestPath() #Reevalute scores with best alignment reset to 0
			print(self)
			#Align in local mode, from the new maximum value (sub-optimal)
			yield from self.__align(self._maxScoreRow, self._maxScoreCol)	
	
	
	def __initialize(self, iniGapPenalty, extGapPenalty):
		#Gap penalties
		self._iniGapPenalty = iniGapPenalty
		self._extGapPenalty = iniGapPenalty if extGapPenalty is None else extGapPenalty
		
		#Fill initial values
		for i in range(1, max(len(self._colSeq), len(self._rowSeq))+1):	
			#Global alignment : first line and colunm have initial scores
			if self._alignMode== "global":
				val = self._iniGapPenalty + self._extGapPenalty*(i-1)
				if i <= len(self._colSeq): #First row
					self._alignMatrix[0][i] = val
					self._rowGapMatrix[0][i] = val 
					self._colGapMatrix[0][i] = val
				
				if i <= len(self._rowSeq): #First column
					self._alignMatrix[i][0] = val
					self._rowGapMatrix[i][0] = val 
					self._colGapMatrix[i][0] = val
					
			#Non-local alignment : first line and column have initial origins
			if self._alignMode != "local":
				if i <= len(self._colSeq): #First row
					self._originMatrix[0][i] = "L"
				
				if i <= len(self._rowSeq): #First column
					self._originMatrix[i][0] = "T"
			
				
		self._alignSeqA = Sequence() #Alignment sequences
		self._alignSeqB = Sequence()
		
		self._maxAlignScore = 0
		self._bestAlignPath = [] #sequence of indexes for alignment with best score
		self._currentAlignPath = []
		
		self._maxScoreRow = 0 #Row of maximum value
		self._maxScoreCol = 0 #Column of maximum value 
		
		#Fill all matrix
		for row in range(1, len(self._rowSeq)+1):
			for col in range(1, len(self._colSeq)+1):
				self.__fill(row, col)
	
	
	def __fill(self, row, col):
		#Row gap matrix
		initialGapScore = self._alignMatrix[row-1][col] + self._iniGapPenalty
		extendedGapScore = self._rowGapMatrix[row-1][col] + self._extGapPenalty
		self._rowGapMatrix[row][col] = max(initialGapScore, extendedGapScore)
		
		#Column gap matrix
		initialGapScore = self._alignMatrix[row][col-1] + self._iniGapPenalty
		extendedGapScore = self._colGapMatrix[row][col-1] + self._extGapPenalty
		self._colGapMatrix[row][col] = max(initialGapScore, extendedGapScore)
		
		#Align matrix
		deleteScore = self._rowGapMatrix[row][col]
		insertScore = self._colGapMatrix[row][col]
		matchScore = self._alignMatrix[row-1][col-1] + \
			self._scoreMatrix.getScore(self._rowSeq[row-1], self._colSeq[col-1])
			
		maxScore = 0
		if self._alignMode == "local": #Local alignment
			self._alignMatrix[row][col] = maxScore = max(insertScore, deleteScore, matchScore, 0)
		else:	#Global alignment
			self._alignMatrix[row][col] = maxScore = max(insertScore, deleteScore, matchScore)
			
		#Max align values
		if maxScore > self._alignMatrix[self._maxScoreRow][self._maxScoreCol]:
			self._maxScoreRow = row
			self._maxScoreCol = col
			
		#Origin matrix
		if insertScore == maxScore:
			self._originMatrix[row][col] += "L"
		if matchScore == maxScore:
			self._originMatrix[row][col] += "D"
		if deleteScore == maxScore:
			self._originMatrix[row][col] += "T"
	
	
	def __clearBestPath(self):
		self._maxAlignScore = 0
		self._maxScoreRow = 0 #Row of maximum value
		self._maxScoreCol = 0 #Column of maximum value
		
		for row, col in self._bestAlignPath[::-1]:
			self._alignMatrix[row][col] = 0
			self._colGapMatrix[row][col] = 0
			self._rowGapMatrix[row][col] = 0
			
			for furtherRow in range(row+1, len(self._rowSeq)+1):
				if (furtherRow, col) not in self._bestAlignPath:
					self._originMatrix[furtherRow][col] = ""
					self.__fill(furtherRow, col)
			
			for furtherCol in range(col+1, len(self._colSeq)+1):
				if (row, furtherCol) not in self._bestAlignPath:
					self._originMatrix[row][furtherCol] = ""
					self.__fill(row, furtherCol)
			
		
	
	def __align(self, i, j, score=0):
		#Local alignments are complete when we reach a null score
		#Global alignments are complete when we reach the beginning of the matrix
		if (self._alignMode == "local" and self._alignMatrix[i][j]==0) \
		or (self._alignMode != "local" and i==0 and j==0):
			yield Sequence(self._alignSeqA), Sequence(self._alignSeqB), score #yield copies of aligned sequences and score
			
			if score > self._maxAlignScore: #Best alignment is remembered for subobtimal lookup
				self._maxAlignScore = score
				self._bestAlignPath = deepcopy(self._currentAlignPath)
			
		else:
			score += self._alignMatrix[i][j]
			for origin in self._originMatrix[i][j]:
				self._currentAlignPath.append((i,j))
				
				if origin == "T": #top
					self._alignSeqA.insert("gap", 0)
					self._alignSeqB.insert(self._rowSeq[i-1], 0)
					yield from self.__align(i-1, j, score)
					
				elif origin == "D": #diagonal
					self._alignSeqA.insert(self._colSeq[j-1], 0)
					self._alignSeqB.insert(self._rowSeq[i-1], 0)
					yield from self.__align(i-1, j-1, score)
				
				elif origin == "L": #left
					self._alignSeqA.insert(self._colSeq[j-1], 0)
					self._alignSeqB.insert("gap", 0)
					yield from self.__align(i, j-1, score)
				else:
					raise ValueError("Origin must be T (top), D (diagonal) or L (left)")
				
				self._currentAlignPath.pop()
				self._alignSeqA.delete(0)
				self._alignSeqB.delete(0)

		
	#Representation
	def __repr__(self):
		sepSize = 5
		rep = {"":"", "L": "_   ", "D": " \\  ", "T": "   |", "LD": "_\\  ","DT": " \\ |", "LT": "_  |", "LDT": "_\\ |"}
		#Amino Acid column sequence
		print(" "*2*sepSize, end="")
		for aa in self._colSeq :
			print('{a!s:<{w}}'.format(a=aa, w=sepSize), end="")
		
		#First line of values (no amino acid)
		print("\n", end=" ")
		for origin in self._originMatrix[0]:
			print('{o:>{w}}'.format(o=rep[origin], w=sepSize), end="")			
		print("\n", end=" "*sepSize)
		for value in self._alignMatrix[0]:
			print('{v:<{w}}'.format(v=value, w=sepSize), end="")
		
		#Amino Acid line sequence and values
		for values, origins, aa in zip(self._alignMatrix[1:], self._originMatrix[1:], self._rowSeq):
			print("\n", end=" ")
			for origin in origins:
				print('{o:>{w}}'.format(o=rep[origin], w=sepSize), end="")
			print()
			print('{a!s:<{w}}'.format(a=aa, w=sepSize), end="")
			for value in values:
				print('{v:<{w}}'.format(v=value, w=sepSize), end="")
		return ""
		

s = ScoreMatrix(r"..\Resources\blosum\blosum62.iij")
print(s)

a = AlignMatrix(s, Sequence("ISALIGNED"), Sequence("THISLINE"))
for seqA, seqB, score in a.globalAlign(-8, -8, True):
	print("A : ", seqA)
	print("B : ", seqB)
	print("Score : ", score)


"""
seq = Sequence.loadFasta(r"..\Resources\fasta\maguk-sequences.fasta")
		
print(seq[0])
print(seq[1])
a = AlignMatrix(s, seq[0], seq[1])

res = []
for seqA, seqB in a.globalAlign(4, 1):
	print(seqA[0:20], seqB[0:20], sep="\n")
	break
"""


		