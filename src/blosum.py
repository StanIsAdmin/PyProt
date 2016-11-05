from sequence import Sequence

def createFromFasta(path):
	sequences = Sequence.loadFasta(path)
	for seq in sequences:
		print(seq.getDescription())
		
createFromFasta(r"C:\Users\mytra\Documents\GitHub\BioInfo\Resources\fasta\PDZ-sequences.fasta")