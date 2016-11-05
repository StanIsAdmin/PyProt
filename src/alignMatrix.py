from scoreMatrix import ScoreMatrix
from sequence import Sequence
from aminoAcid import AminoAcid
from copy import deepcopy

class AlignedSequences:
	def __init__(self, seqA, seqB, score, identity, gaps, alignmentType):
		self.seqA = seqA
		self.seqB = seqB
		self.score = score
		self.identity = identity
		self.identityPercent = 100*(identity/len(seqA))
		self.gaps = gaps
		self.alignmentType = alignmentType
		self.chunkSize = 80
		
	def __repr__(self):
		res = []
		res.append("========== Alignment ==========")
		res.append("Type     : " + self.alignmentType)
		res.append("Score    : " + str(self.score))
		res.append("Identity : " + str(self.identity) + " ({0:.2f}%)".format(self.identityPercent))
		res.append("Gaps     : " + (str(self.gaps) if self.gaps > 0 else "None"))
		res.append("")
		
		if "" not in (self.seqA.getDescription(), self.seqB.getDescription()):
			res.append("Upper Sequence : " + self.seqA.getDescription())
			res.append("Lower Sequence : " + self.seqB.getDescription())
			res.append("")
		
		listA, listB, listSep = [], [], []
		for a, b in zip(self.seqA, self.seqB):
			listA.append(str(a))
			listB.append(str(b))
			if a==b:
				listSep.append("|")
			elif a.isGap() or b.isGap():
				listSep.append(" ")
			else:
				listSep.append(":")
			
			
			if len(listA) == self.chunkSize:
				res.append("".join(listA))
				res.append("".join(listSep))
				res.append("".join(listB))
				res.append("\n")
				listA, listB, listSep = [], [], []
		
		if len(listA) > 0:
				res.append("".join(listA))
				res.append("".join(listSep))
				res.append("".join(listB))
		res.append("")

		return "\n".join(res)

		

class AlignMatrix:
	"""
	Represents an alignment matrix, used to determine the alignment score between two sequences
	"""
	
	def __init__(self, scoreMatrix):
		"""
		Creates an AlignMatrix object that uses scoreMatrix as a scoring system between amino acids.
		"""
		self._scoreMatrix = scoreMatrix #scoring matrix to use
		
		self._colSeq = None #copies of sequences to align
		self._rowSeq = None
		
		self._alignMatrix = None #Used for alignment scores between sequences
		self._rowGapMatrix = None #used for affine gap penalty
		self._colGapMatrix = None
		self._originMatrix = None #Used for backtracking purposes
		
		self._alignMode = "" #Align mode (global, semiglobal, local)
		
		self._maxAlignScore = 0
		self._maxScoreRow = 0 #Index of maximum value from _alignMatrix
		self._maxScoreCol = 0
		self._currentAlignPath = []
		self._bestAlignPath = [] #sequence of origins for alignment with best score
		
		self._iniGapPenalty = 0 #Initial gap penalty
		self._extGapPenalty = 0 #Extended gap penalty
		
		self._alignedSeqA = None #Aligned sequences (result of alignment)
		self._alignedSeqB = None
		
	def firstGlobalAlign(self, seqA, seqB, iniGapPenalty, extGapPenalty=None):
		for result in self.globalAlign(seqA, seqB, iniGapPenalty, extGapPenalty):
			return result
	
	def globalAlign(self, seqA, seqB, iniGapPenalty, extGapPenalty=None, semiGlobal=False):
		if semiGlobal:
			self._alignMode = "semiglobal"
		else:
			self._alignMode = "global"
		
		#Initialize all data structures
		self.__initialize(seqA, seqB, iniGapPenalty, extGapPenalty)
		
		#Align in global mode, from the bottom right
		yield from self.__align(len(self._rowSeq), len(self._colSeq))
	
	
	def firstLocalAlign(self, seqA, seqB, iniGapPenalty, extGapPenalty=None):
		for result in self.localAlign(seqA, seqB, iniGapPenalty, extGapPenalty):
			return result
	
	def localAlign(self, seqA, seqB, iniGapPenalty, extGapPenalty=None, subOptimal=False):
		self._alignMode = "local"
		
		#Initialize all data structures
		self.__initialize(seqA, seqB, iniGapPenalty, extGapPenalty)
		
		#Align in local mode, from the maximum value (optimal)
		yield from self.__align(self._maxScoreRow, self._maxScoreCol)
		
		#If suboptimal alignments are requested
		if subOptimal:
			self.__clearBestPath() #Reevalute scores with best alignment reset to 0
			#Align in local mode, from the new maximum value (sub-optimal)
			yield from self.__align(self._maxScoreRow, self._maxScoreCol)	
	
	
	def __initialize(self, seqA, seqB, iniGapPenalty, extGapPenalty):
		#Sequences
		if len(seqA)==0 or len(seqB)==0:
				raise ValueError("seqA and seqB cannot be empty")
		self._colSeq = Sequence(seqA)
		self._rowSeq = Sequence(seqB)
		
		#Alignments
		self._alignSeqA = Sequence()
		self._alignSeqB = Sequence()
		
		self._maxAlignScore = 0 #Maximum score found while aligning
		self._maxScoreRow = 0 #Row of maximum score
		self._maxScoreCol = 0 #Column of maximum score 
		
		self._currentAlignPath = [] #indexes of the current alignment
		self._bestAlignPath = [] #indexes of the alignment with the best score
		
		#Gap penalties
		self._iniGapPenalty = iniGapPenalty
		self._extGapPenalty = iniGapPenalty if extGapPenalty is None else extGapPenalty
		
		#Matrices
		self._alignMatrix = [[0 for i in range(len(self._colSeq)+1)] for j in range(len(self._rowSeq)+1)]
		self._rowGapMatrix = deepcopy(self._alignMatrix)
		self._colGapMatrix = deepcopy(self._alignMatrix)
		self._originMatrix = [["" for i in range(len(self._colSeq)+1)] for j in range(len(self._rowSeq)+1)]
		
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
			
		
		#Fill all matrices
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
			
		
	
	def __align(self, i, j, score=0, identity=0, gaps=0):
		#Local alignments are complete when we reach a null score
		#Global alignments are complete when we reach the beginning of the matrix
		if (self._alignMode == "local" and self._alignMatrix[i][j]==0) \
		or (self._alignMode != "local" and i==0 and j==0):
			#Create AlignedSequences object as result
			result = AlignedSequences(Sequence(self._alignSeqA), Sequence(self._alignSeqB),\
			score, identity, gaps, self._alignMode)
			
			#Remember best alignment for subobtimal lookup
			if score > self._maxAlignScore: 
				self._maxAlignScore = score
				self._bestAlignPath = deepcopy(self._currentAlignPath)
			
			yield result
			
		else:
			score += self._alignMatrix[i][j] #TODO: check if score is sum or just last value
			for origin in self._originMatrix[i][j]:
				self._currentAlignPath.append((i,j))
				
				if origin == "T": #top
					gaps += 1
					self._alignSeqA.insert("gap", 0)
					self._alignSeqB.insert(self._rowSeq[i-1], 0)
					yield from self.__align(i-1, j, score, identity, gaps)
					
				elif origin == "D": #diagonal
					if self._colSeq[j-1] == self._rowSeq[i-1]:
						identity += 1
					self._alignSeqA.insert(self._colSeq[j-1], 0)
					self._alignSeqB.insert(self._rowSeq[i-1], 0)
					yield from self.__align(i-1, j-1, score, identity, gaps)
				
				elif origin == "L": #left
					gaps += 1
					self._alignSeqA.insert(self._colSeq[j-1], 0)
					self._alignSeqB.insert("gap", 0)
					yield from self.__align(i, j-1, score, identity, gaps)
				
				else:
					raise ValueError("Origin must be T (top), D (diagonal) or L (left)")
				
				self._currentAlignPath.pop()
				self._alignSeqA.delete(0)
				self._alignSeqB.delete(0)

		
	#Representation
	def __repr__(self):
		result = []
		sepSize = 5
		rep = {"":"", "L": "_   ", "D": " \\  ", "T": "   |", "LD": "_\\  ","DT": " \\ |", "LT": "_  |", "LDT": "_\\ |"}
		#Amino Acid column sequence
		result.append(" "*2*sepSize)
		for aa in self._colSeq :
			result[-1] += '{a!s:<{w}}'.format(a=aa, w=sepSize)
		
		#First line of values (no amino acid)
		result.append(" ")
		for origin in self._originMatrix[0]:
			result[-1] += '{o:>{w}}'.format(o=rep[origin], w=sepSize)
		
		result.append(" "*sepSize)
		for value in self._alignMatrix[0]:
			result[-1] += '{v:<{w}}'.format(v=value, w=sepSize)
		
		#Amino Acid line sequence and values
		for values, origins, aa in zip(self._alignMatrix[1:], self._originMatrix[1:], self._rowSeq):
			result.append(" ")
			for origin in origins:
				result[-1] += '{o:>{w}}'.format(o=rep[origin], w=sepSize)
			result.append("")
			result[-1] += '{a!s:<{w}}'.format(a=aa, w=sepSize)
			for value in values:
				result[-1] += '{v:<{w}}'.format(v=value, w=sepSize)
		
		return "\n".join(result)
		