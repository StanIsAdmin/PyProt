from math import sqrt, log, ceil

from pyprot.protein import AminoAcid, Sequence, loadFasta


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
    Creates a profile for a series of aligned sequences, and gives a score to each AA subsitution in a given column.
    """

    def __init__(self, description=""):
        self.description = description
        self.seqCount = 0  # total number of sequences
        self.size = None  # all sequences have the same size
        self.aaDistribution = None  # amino acid distribution
        self.aaCount = None
        self.gapPenalties = None

    def add(self, sequence):
        # check sequence size
        if self.size is None:
            self.size = len(sequence)
            self.aaDistribution = [{} for i in range(self.size)]
            self.aaCount = [0 for i in range(self.size)]
            self.gapPenalties = [0 for i in range(self.size + 1)]

        assert (len(sequence) == self.size)

        # update amino acid count for each column
        for index in range(self.size):
            if not sequence[index].isGap():
                self.aaCount[index] += 1
                try:
                    self.aaDistribution[index][sequence[index]] += 1
                except:
                    self.aaDistribution[index][sequence[index]] = 1

        # increase sequence count
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


def belongs(outSequence, group, minMatches):
    """
    Returns True if there's at least 'minMatches' matches between outSequence and any sequence 
    from 'group'; False otherwise.
    Size of all sequences is assumed to be equal.
    """
    for inSequence in group:
        matches = 0
        for i in range(len(outSequence)):
            if outSequence[i] == inSequence[i]:
                matches += 1
            if matches >= minMatches:  # As soon as a match is found
                return True
    return False


def makeGroupsFromFasta(path, requiredIdentityPercent):
    """
    Loads Sequences from file at 'path' and separate them in groups with an identity of 
    at least 'requiredIdentityPercent'.
    Returns a list representing the groups as lists of Sequence objects.
    """
    sequences = [seq for seq in loadFasta(path)]
    groups = [[sequences[0]]]  # First sequence is assigned to first group

    seqSize = len(sequences[0])  # Size of the sequences

    # Number of matches required to achieve requiredIdentityPercent
    minMatches = ceil((requiredIdentityPercent / 100) * seqSize)

    # For each outSequence not yet in a group
    for outSequence in sequences[1:]:
        groupFound = False  # has a group been found ?
        groupIndex = 0  # index of the group we're looking in

        # Look for a group where outSequence belongs
        while (not groupFound) and groupIndex < len(groups):
            if belongs(outSequence, groups[groupIndex], minMatches):
                groups[groupIndex].append(outSequence)
                groupFound = True
            else:
                groupIndex += 1  # Move on to next group

        # If no group works, create a new one
        if not groupFound:
            groups.append([outSequence])

    return groups


def valueDictsFromGroups(groups):
    """
    Transforms each group from 'groups' into a list of dictionaries, one for each Sequence column.
    The dictionaries map each AminoAcid found in that column to their count.
    Returns the list of dictionaries and a list of the size (in Sequences) of their groups.
    """
    seqSize = len(groups[0][0])
    groupCount = len(groups)

    groupValues = []
    groupSizes = []

    for group in groups:  # For each group
        groupAAs = []
        groupSize = len(group)  # Size of group (n°sequences)

        for col in range(seqSize):  # For each column
            groupCol = {}

            for seq in group:  # For each sequence in group
                try:
                    groupCol[seq[col]] += 1  # Increment count
                except:
                    groupCol[seq[col]] = 1

            groupAAs.append(groupCol)
        groupValues.append(groupAAs)
        groupSizes.append(groupSize)
    return groupValues, groupSizes


def getFrequencies(groupValues, groupSizes):
    """
    Evaluates the frequencies of AminoAcids within columns of groups in 'groupValues'.
    Frequencies are weighted according to group sizes in 'groupSizes'.
    Returns two dictionaries and a number:
        -'freqPairs' maps pairs of AminoAcids to their frequencies
        -'freqSingle' maps single AminoAcids to their frequencies
        -'freqSum' is the sum of all frequencies
    """
    seqSize = len(groupValues[0])  # Size of the Sequences
    groupCount = len(groupSizes)  # Number of groups

    freqPairs = {}  # frequencies of amino acid pairs (fAB)

    # Frequencies of single amino acids (fA)
    freqSingle = {AminoAcid(aa): 0 for aa in AminoAcid.getAllNames()}

    freqSum = 0  # Sum of frequencies  sum(fAB)

    for col in range(seqSize):  # Each column
        for groupAIndex in range(groupCount - 1):  # Each groupA
            groupA = groupValues[groupAIndex]
            groupASize = groupSizes[groupAIndex]

            for groupBIndex in range(groupAIndex + 1, groupCount):  # Each further groupB
                groupB = groupValues[groupBIndex]
                groupBSize = groupSizes[groupBIndex]

                for aaA, aaACount in groupA[col].items():  # Each AA from groupA
                    aaAFreq = aaACount / groupASize  # Its frequency within groupA

                    for aaB, aaBCount in groupB[col].items():  # Each AA from groupB
                        aaBFreq = aaBCount / groupBSize  # Its frequency within groupB

                        aaPairFreq = aaAFreq * aaBFreq  # Pair frequency
                        freqSum += aaPairFreq  # Sum of all frequencies
                        freqSingle[aaA] += aaPairFreq / 2
                        freqSingle[aaB] += aaPairFreq / 2

                        # Index is unique to this pair
                        pairIndex = (aaA, aaB) if aaA > aaB else (aaB, aaA)
                        try:
                            freqPairs[pairIndex] += aaPairFreq
                        except:
                            freqPairs[pairIndex] = aaPairFreq

    return freqPairs, freqSingle, freqSum


def sumFrequenciesToProb(freqPairsList, freqSingleList, freqSumList):
    """
    Sums all frequencies in provided lists, and transforms them to probabilities according
    to the sum of frequencies found in 'freqSumList'.
    """
    fSum = sum(freqSumList)  # Absolute sum of all frequencies
    probPairs, probSingle = {}, {}  # Probabilities for pairs and single AAs

    for freqPairs, freqSingle in zip(freqPairsList, freqSingleList):
        for key, value in freqPairs.items():
            # Sum all frequencies for matching AA pairs, divided by fSum
            try:
                probPairs[key] += value / fSum
            except:
                probPairs[key] = value / fSum
        for key, value in freqSingle.items():
            # Sum all frequencies for matching AAs, divided by fSum
            try:
                probSingle[key] += value / fSum
            except:
                probSingle[key] = value / fSum

    return probPairs, probSingle


def blosumFromProbabilities(probPairs, probSingle, requiredIdentityPercent):
    """
    Fills and returns a Score according to the BLOSUM algorithm, from the probabilities
    of AA pairs and singletons provided in 'probPairs' and 'probSingle'.
    """
    # Create empty Score, ignoring AAs B,Z,J,U,O
    scoreMatrix = ScoreMatrix("", "BLOSUM{}".format(requiredIdentityPercent), "BZJUO")

    for key, qAB in probPairs.items():
        # qAB is the evolutionary probability of the AA pair (A,B)
        aaA, aaB = key  # Both amino acids
        pA, pB = probSingle[aaA], probSingle[aaB]  # Their single probabilities

        # eAB is the random probability of the AA pair (A, B) given their single probabilities
        eAB = pA * pA if aaA == aaB else 2 * pA * pB

        # The BLOSUM score for this pair is the log-odds-ratio of evolutionary and random prob.
        sAB = int(round(2 * log(qAB / eAB, 2)))
        scoreMatrix.setScore(aaA, aaB, sAB)  # Fill the matrix

    return scoreMatrix


def blosumFromFasta(requiredIdentityPercent, *filepaths):
    """
    Creates and returns a Score for all sequences in the provided 'filepaths',
    using the BLOSUM approach with an identity of at least 'requiredIdentityPercent'.
    Each file is grouped independently and only then their weighted probabilities are merged.
    """
    # Results for each different files are first stored in lists
    freqPairsList, freqSingleList, freqSumList = [], [], []

    for path in filepaths:  # for each .fasta file
        groups = makeGroupsFromFasta(path, requiredIdentityPercent)  # groups
        groupValues, groupSizes = valueDictsFromGroups(groups)  # groups as dicts
        freqPairs, freqSingle, freqSum = getFrequencies(groupValues, groupSizes)  # frequencies

        freqPairsList.append(freqPairs)  # Append results
        freqSingleList.append(freqSingle)
        freqSumList.append(freqSum)

    # Merge (sum) results together
    probPairs, probSingle = sumFrequenciesToProb(freqPairsList, freqSingleList, freqSumList)
    # Create the BLOSUM matrix
    blosum = blosumFromProbabilities(probPairs, probSingle, requiredIdentityPercent)
    print(blosum)
