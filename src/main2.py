from dssp import DSSP

d = DSSP(r"../resources/dataset/dssp/1A2P.dssp")
seq, struct = d.getSequenceStructure("A")
print(">", d.identifier, d.protein, d.organism)
print(seq)
print(struct)