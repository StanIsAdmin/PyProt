from sequence import Sequence, loadFasta
from score import PSSM
from align import Align

pssm = PSSM("WW domain")
for seq in loadFasta(r"../resources/fasta/msaresults-MUSCLE.fasta"):
	pssm.add(seq)
pssm.setGapPenalty(4)

al = Align(pssm)

for toalign in loadFasta(r"../resources/fasta/test.fasta"):
	for aligned in al.multiAlign(toalign):
		print(aligned)