from copy import deepcopy
from math import sqrt
from aminoacid import AminoAcid
from sequence import Sequence, loadFasta
from score import Score


# AA frequencies for complete UniProt database
# from http://web.expasy.org/docs/relnotes/relstat.html, "AMINO ACID COMPOSITION"
uniprob = {
	AminoAcid("Ala") : .0826,
	AminoAcid("Gln") : .0393,
	AminoAcid("Leu") : .0965,
	AminoAcid("Ser") : .0660,
	AminoAcid("Arg") : .0553,
	AminoAcid("Glu") : .0674,
	AminoAcid("Lys") : .0582,
	AminoAcid("Thr") : .0535,
	AminoAcid("Asn") : .0406,
	AminoAcid("Gly") : .0708,
	AminoAcid("Met") : .0241,
	AminoAcid("Trp") : .0109,
	AminoAcid("Asp") : .0546,
	AminoAcid("His") : .0227,
	AminoAcid("Phe") : .0386,
	AminoAcid("Tyr") : .0292,
	AminoAcid("Cys") : .0137,
	AminoAcid("Ile") : .0593,
	AminoAcid("Pro") : .0472,
	AminoAcid("Val") : .0687,
	
}	


class MAlign:
	"""
	'Multiple Alignments'
	Aligns multiple sequences together.
	"""
	def __init__(self, scoreMatrix, sequenceLength):
		self.scoreMatrix = scoreMatrix #scoring matrix used for alignment
		self.seqLen = sequenceLength #all sequences have the same size
		self.aaDistribution = [{} for i in range(self.seqLen)] #amino acid distribution
		
		self.sequences = [] #aligned sequences
		self.seqCount = 0 #total number of sequences
		
	
	def align(self, sequence):
		#pseudocounts
		alpha = self.seqCount - 1
		beta = sqrt(self.seqCount)
		alphaplusbeta = alpha + beta
		
		#scoring system
		scoreMatrix = [[0 for i in range(self.seqLen)] for j in range(self.seqLen)]
		
		index = -1
		for aa in sequence:
			index += 1
			
			#random probability of amino acid
			try:
				p_aa = uniprob[aa]
			except:
				p_aa = 0.001
			
			#evolutionary probability of amino acid
			try:
				f_aa = self.aaDistribution[index][aa] / self.seqCount
			except:
				f_aa = 0
				
			q_aa = (alpha * f_aa + beta * p_aa) / alphaplusbeta
			m_aa = q_aa / p_aa
			
	
	def __add(self, sequence):
		#check sequence size
		assert(len(sequence) == self.seqLen)
			
		#update amino acid count for each column
		for index in range(self.seqLen):
			try:
				self.aaDistribution[index][sequence[index]] += 1
			except:
				self.aaDistribution[index][sequence[index]] = 1
		
		#add sequence
		self.sequences.append(sequence)
		self.seqCount += 1
