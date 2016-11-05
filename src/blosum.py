from sequence import Sequence
from scoreMatrix import ScoreMatrix

def createFromFasta(path):
	sequences = Sequence.loadFasta(path)
	
	
	
	scoreMat = ScoreMatrix()
		
createFromFasta(r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\PDZ-sequences.fasta")