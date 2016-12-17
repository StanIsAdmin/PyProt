from dssp import DSSP
from structure import GOR3

"""
with open(r"../resources/dataset/CATH_info.txt") as infoFile:
	with open(r"../resources/dataset/CATH_info-PARSED.txt", 'w') as outFile:
		for line in infoFile.readlines():
			d = DSSP(r"../resources/dataset/dssp/" + line[0:4] + ".dssp")
			description = "> " + d.identifier + "|" + d.protein + "|" + d.organism
			seq, struct = d.getSequenceStructure(line[4])
			
			outFile.writelines(l + "\n" for l in [description,seq,struct])
		
with open(r"../resources/dataset/CATH_info_test.txt") as infoFile:
	with open(r"../resources/dataset/CATH_info_test-PARSED.txt", 'w') as outFile:
		for line in infoFile.readlines():
			d = DSSP(r"../resources/dataset/dssp_test/" + line[0:4] + ".dssp")
			description = "> " + d.identifier + "|" + d.protein + "|" + d.organism
			seq, struct = d.getSequenceStructure(line[4])
			
			outFile.writelines(l + "\n" for l in [description,seq,struct])
"""

p = GOR3()

with open(r"../resources/dataset/CATH_info-PARSED.txt") as inFile:
	index = 0
	sequence = ""
	for line in inFile.readlines():
		line = line.strip().upper()
		if not (line=="" or line[0]==">"):
			#Line is a sequence
			if index % 2 == 0:
				sequence = line
			#Line is a structure
			else:
				p.train(sequence, line)
				sequence = ""
			index += 1
			
s = p.predict("RPYACPVESCDRRFSRSADLTRHIRIHTGQKPFQCRICMRNFSRSDHLTTHIRTHTGEKPFACDICGRKFARSDERKRHTKIHLR")
realS = "CCEECCCTTCCCEECCHHHHHHHTHHHHTCCCEECTTTCCEECCHHHHHHHHHHHHCCCCEECTTTCCEECCHHHHHHHHHHHCC"
print(s)
print(realS)
result = ""
accuracy = 0
for s1,s2 in zip(s, realS):
	if s1==s2:
		result += ":"
		accuracy += 1
	else:
		result += " "
print(result)
print(accuracy / len(s))
		