

class DSSP:
    def __init__(self, filePath):

        # Interesting columns : (start index, end index)
        self.columns = ("RESIDUE", "AA", "STRUCTURE")
        self.residues = []

        # Metadata
        self.identifier = ""
        self.protein = ""
        self.organism = ""

        # Parsing
        with open(filePath, 'r') as dsspFile:
            columnIndex = {col: (0, 0) for col in self.columns}
            lineIsData = False

            for line in dsspFile.readlines():

                if lineIsData:
                    data = []
                    for column in self.columns:
                        start, end = columnIndex[column]
                        data.append(line[start:end])

                    self.residues.append(data)

                else:
                    if line.strip()[0] == "#":
                        lineIsData = True
                        for column in self.columns:
                            startIndex = line.find(column)
                            endIndex = startIndex + len(column)
                            endIndex = endIndex + (len(line[endIndex:]) - len(line[endIndex:].lstrip())) - 1
                            columnIndex[column] = (startIndex, endIndex)
                    elif line.startswith("HEADER"):
                        self.identifier = line.split()[-2]

                    elif line.startswith("COMPND"):
                        self.protein = line.split(":")[1].split(";")[0].strip()
                        self.protein = " ".join(self.protein.split())

                    elif line.startswith("SOURCE"):
                        self.organism = line.split(":")[1].split(";")[0].strip()
                        self.organism = " ".join(self.organism.split())

    def __repr__(self):
        res = []
        for values in self.residues:
            res.append(str(values))
        return "\n".join(res)

    def getSequenceStructure(self, chain):
        structs = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E", "T": "T", "C": "C", "S": "C", " ": "C"}
        sequence = []
        structure = []

        for residue in self.residues:
            if residue[0][-1] == chain:
                sequence.append(residue[1][0])
                structure.append(structs[residue[2][0]])

        return "".join(sequence), "".join(structure)
