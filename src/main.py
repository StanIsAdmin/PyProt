from aminoacid import AminoAcid
from sequence import Sequence
from score import Score
from talign import TAlign, TAligned

a = AminoAcid("A")
s = Sequence([a, a, "B"])
score = Score()
talign = TAlign(score)