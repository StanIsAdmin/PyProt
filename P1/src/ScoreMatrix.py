class ScoreMatrix:
    """
    Represents a scoring matrix, used to determine the score between two Amino Acids
    """
    
    def __init__(self, path=""):
        self._matrix = []
        lineSize = 1
        
        if path != "":
            with open(path, 'r') as file:
                for line in file:
                    if line[0] != "#" and line.strip()[0] != "A":
                        self._matrix.append([int(v) for v in line.split()])
        else:
            for aa in AminoAcid.getAllNames():
                self._matrix.append([0 for i in range(lineSize)])
                lineSize += 1
            
    def setScore(self, score, aa1, aa2):
        id1 = aa1.getId()
        id2 = aa2.getId()
        if id1 > id2:
            self._matrix[id1][id2] = score
        else:
            self._matrix[id2][id1] = score
        
    def getScore(self, aa1, aa2):
        id1 = aa1.getId()
        id2 = aa2.getId()
        if id1 > id2:
            return self._matrix[id1][id2]
        else:
            return self._matrix[id2][id1]

                
        
            
        
ScoreMatrix(r"Resources\blosum\blosum30.iij")
