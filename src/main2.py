from dssp import DSSP


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
