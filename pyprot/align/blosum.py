from math import log, ceil

from pyprot.base.aminoacid import AminoAcid
from pyprot.data.fasta import getSequencesFromFasta
from pyprot.align.score import ScoreMatrix


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
    Returns a list representing the groups as lists of Protein objects.
    """
    sequences = [seq for seq in getSequencesFromFasta(path)]
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
    Transforms each group from 'groups' into a list of dictionaries, one for each Protein column.
    The dictionaries map each AminoAcid found in that column to their count.
    Returns the list of dictionaries and a list of the size (in Sequences) of their groups.
    """
    seqSize = len(groups[0][0])
    groupCount = len(groups)

    groupValues = []
    groupSizes = []

    for group in groups:  # For each group
        groupAAs = []
        groupSize = len(group)  # Size of group (n sequences)

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
    Creates and returns a ScoreMatrix for all sequences in the provided 'filepaths',
    using the BLOSUM approach with an identity of at least 'requiredIdentityPercent'.
    Each file is grouped independently and only then their weighted probabilities are merged.
    :param requiredIdentityPercent: 
    :param filepaths: 
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
