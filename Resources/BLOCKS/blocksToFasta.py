
blocksfile = r"PDZB.blocks.txt"
fastafile = r"PDZB.fasta"
with open(blocksfile,'r') as f:
	with open(fastafile, 'w') as c:
		for line in f:
			line = line.strip()
			if line != "":
				name, data = line.split(")")
				data = data.split()[0].strip()
				name1, name2 = name.split("(")
				newname = ">" + name1.strip() + "|" + name2.strip() + "-" + str(int(name2)+len(data)-1) + "\n"
				c.write(newname)
				c.write(data + "\n")