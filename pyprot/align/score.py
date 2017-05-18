from math import sqrt, log

from pyprot.base.aminoacid import AminoAcid
from pyprot.base.sequence import Sequence


class ScoreMatrix:
    """
    Represents a scoring matrix, used to determine the score between two Amino Acids
    """

    def __init__(self, path="", description="", ignore=None):
        """
        Creates a Score object.
        If 'path' is provided, loads the Score values from an iij file.
        Otherwise, creates a Score for all possible AminoAcids with values 0.
        """
        self._description = description
        self._ignore = Sequence(ignore)
        self._matrix = []
        self._aaOrder = {}
        self._aaSequence = Sequence()

        # If path is provided, load directly from iij file
        if path != "":
            with open(path, 'r') as file:
                foundAAOrder = False  # Have we found the line with the amino acid values and order yet?
                for line in file:
                    if line[0] != "#":  # Comments

                        if not foundAAOrder:  # Read aa values and order
                            for aa in line.split():
                                self._aaSequence.extend(aa)
                            self._aaOrder = {aa: index for aa, index in
                                             zip(self._aaSequence, range(len(self._aaSequence)))}
                            foundAAOrder = True
                        else:  # Read matrix values
                            self._matrix.append([int(v) for v in line.split()])

        # Otherwise initialize matrix with 0
        else:
            lineSize = 1
            for aa in AminoAcid.getAllNames():
                if AminoAcid(aa) not in self._ignore:
                    self._aaSequence.extend(aa)
                    self._aaOrder[self._aaSequence[-1]] = lineSize - 1
                    self._matrix.append([0 for i in range(lineSize)])
                    lineSize += 1

    # Representation
    def __repr__(self):
        """
        Representation.
        """
        sepSize = 4
        result = ["---------- " + self._description + " ----------"]
        for values, aa in zip(self._matrix, self._aaSequence):
            tempstr = '{a!s:<{w}}'.format(a=aa, w=sepSize)
            for value in values:
                tempstr += '{v:<{w}}'.format(v=value, w=sepSize)
            result.append(tempstr)
        tempstr = " " * sepSize
        for aa in self._aaSequence:
            tempstr += '{a!s:<{w}}'.format(a=aa, w=sepSize)
        result.append("")
        result.append(tempstr)
        return "\n".join(result)

    # Scoring
    def setScore(self, aa1, aa2, score):
        """
        Set the score assigned to AminoAcids 'aa1', 'aa2'.
        """
        id1 = self._aaOrder[aa1]
        id2 = self._aaOrder[aa2]
        if id1 > id2:
            self._matrix[id1][id2] = score
        else:
            self._matrix[id2][id1] = score

    def getScore(self, aa1, aa2):
        """
        Get the score assigned to AminoAcids 'aa1', 'aa2'.
        """
        id1 = self._aaOrder[aa1]
        id2 = self._aaOrder[aa2]
        if id1 > id2:
            return self._matrix[id1][id2]
        else:
            return self._matrix[id2][id1]


# AA frequencies for complete UniProt database
# from http://web.expasy.org/docs/relnotes/relstat.html, "AMINO ACID COMPOSITION"
uniprob = {
    AminoAcid("Ala"): .0826,
    AminoAcid("Gln"): .0393,
    AminoAcid("Leu"): .0965,
    AminoAcid("Ser"): .0660,
    AminoAcid("Arg"): .0553,
    AminoAcid("Glu"): .0674,
    AminoAcid("Lys"): .0582,
    AminoAcid("Thr"): .0535,
    AminoAcid("Asn"): .0406,
    AminoAcid("Gly"): .0708,
    AminoAcid("Met"): .0241,
    AminoAcid("Trp"): .0109,
    AminoAcid("Asp"): .0546,
    AminoAcid("His"): .0227,
    AminoAcid("Phe"): .0386,
    AminoAcid("Tyr"): .0292,
    AminoAcid("Cys"): .0137,
    AminoAcid("Ile"): .0593,
    AminoAcid("Pro"): .0472,
    AminoAcid("Val"): .0687,

}


class PSSM:
    """
    Position Specific Score Matrix.
    Creates a profile for a series of aligned Sequences, and gives a score to each AA subsitution in a given column.
    """

    def __init__(self, description=""):
        self.description = description
        self.seqCount = 0  # total number of Sequences
        self.size = None  # all Sequences have the same size
        self.aaDistribution = None  # amino acid distribution
        self.aaCount = None
        self.gapPenalties = None

    def add(self, Sequence):
        # check Sequence size
        if self.size is None:
            self.size = len(Sequence)
            self.aaDistribution = [{} for i in range(self.size)]
            self.aaCount = [0 for i in range(self.size)]
            self.gapPenalties = [0 for i in range(self.size + 1)]

        assert (len(Sequence) == self.size)

        # update amino acid count for each column
        for index in range(self.size):
            if not Sequence[index].isGap():
                self.aaCount[index] += 1
                try:
                    self.aaDistribution[index][Sequence[index]] += 1
                except:
                    self.aaDistribution[index][Sequence[index]] = 1

        # increase Sequence count
        self.seqCount += 1

    def getDescription(self):
        return self.description

    def getScore(self, aminoAcid, columnIndex):
        # pseudocounts
        alpha = self.aaCount[columnIndex] - 1
        beta = sqrt(self.seqCount)
        alphaplusbeta = alpha + beta

        # random probability of amino acid
        try:
            p_aa = uniprob[aminoAcid]
        except:
            p_aa = 0.001

        # evolutionary probability of amino acid
        try:
            f_aa = self.aaDistribution[columnIndex][aminoAcid] / self.seqCount
        except:
            f_aa = 0

        q_aa = (alpha * f_aa + beta * p_aa) / alphaplusbeta

        return log(q_aa / p_aa)

    def getGapPenalty(self, columnIndex):
        return self.gapPenalties[columnIndex]

    def setGapPenalty(self, penalty, columnIndex=None):
        if columnIndex is None:
            for i in range(self.size):
                self.gapPenalties[i] = penalty
        else:
            self.gapPenalties[columnIndex] = penalty

    def __len__(self):
        return self.size

    def __repr__(self):
        for i in range(self.size):
            for key, score in self.aaDistribution[i].items():
                print(key, ": ", score, "(", self.getScore(key, i), ")", sep="", end=", ")
            print()
