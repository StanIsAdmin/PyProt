from aminoacid import AminoAcid
from sequence import Sequence, loadFasta
from score import ScoreMatrix, PSSM
from align import Align, Aligned

score = ScoreMatrix(r"../resources/blosum/blosum62.iij")
pssm = PSSM()
for seq in loadFasta(r"../resources/fasta/msaresults-MUSCLE.fasta"):
	pssm.add(seq)
pssm.setGapPenalty(10)

al = Align(pssm)

for toalign in loadFasta(r"../resources/fasta/test.fasta"):
	for aligned in al.multiAlign(toalign):
		#print(aligned.getBottomSequence())
		print(aligned)