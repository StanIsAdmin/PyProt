from aminoacid import AminoAcid
from sequence import Sequence
from score import Score
from talign import TAlign, TAligned
from malign import MAlign

score = Score(r"../resources/blosum/blosum62.iij")
malign = MAlign(score, sequenceLength=3)
malign.align(Sequence("ABC"))
