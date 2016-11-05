from r"C:\Users\mytra\Documents\GitHub\BioInfo\P1\src\sequence.py" import Sequence

def createFromFasta(path):
	sequences = Sequence.loadFasta(path)
	for seq in sequences:
		print(seq.getDescription())