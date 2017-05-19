from copy import deepcopy

from pyprot.base.sequence import Sequence


class Aligned:
    """
    Represents two aligned Sequences with some metadata about the alignment.
    Alignments can be two-sequences alignments or multiple sequence alignments.
    """

    def __init__(self, seqA, seqB, seqAStart, seqBStart, alignType, alignScore, scoreMatrix, isMultiple=False):
        assert (len(seqA) == len(seqB))

        # seqA is always a Sequence, seqB can be a list of ints (if MSA)
        self.seqA = seqA if isinstance(seqA, Sequence) else seqB
        self.seqB = seqB if isinstance(seqA, Sequence) else seqA
        self.seqInter = []  # Sequence separators (.: )

        self.seqAStart = seqAStart if isinstance(seqA, Sequence) else seqBStart
        self.seqBStart = seqBStart if isinstance(seqA, Sequence) else seqAStart

        self.size = len(self.seqA)

        self.alignType = alignType  # Options used for alignemnt
        self.alignScore = round(alignScore, 2)  # Score of alignment

        self.isMultiple = isMultiple  # Is this a M.S.A ?
        self.pssmDescription = None
        if self.isMultiple:
            self.pssmDescription = scoreMatrix.getDescription()

        self.condensed = False  # When true, aligned Sequences are not displayed
        self.chunkSize = 80  # Used for display
        self.sizeIndic = 20  # Same

        self.identity = 0  # Identity (same AAs)
        self.gaps = 0  # Gap count (either AA is a gap)
        self.seqAGaps = 0
        self.seqBGaps = 0
        self.similarity = 0  # Similarity (positive score between AAs)

        if self.isMultiple:
            for aa in self.seqA:
                if aa.isGap():
                    self.gaps += 1
                    self.seqAGaps += 1
        else:
            for aa1, aa2 in zip(self.seqA, self.seqB):
                if aa1.isGap() or aa2.isGap():
                    if aa1.isGap():
                        self.seqAGaps += 1
                    else:
                        self.seqBGaps += 1
                    self.seqInter.append(" ")
                    self.gaps += 1
                elif aa1 == aa2:
                    self.seqInter.append(":")
                    self.identity += 1
                elif scoreMatrix.getScore(aa1, aa2) >= 0:
                    self.seqInter.append(".")
                    self.similarity += 1
                else:
                    self.seqInter.append(" ")
        self.seqAEnd = self.seqAStart + len(self.seqA) - self.seqAGaps
        self.seqBEnd = self.seqBStart + len(self.seqB) - self.seqBGaps

    def getBottomSequence(self):
        return self.seqB

    def __repr__(self):
        """
        Object representation.
        """
        res = []
        if self.isMultiple:
            res.append("---------- Multi-Seq. Alignment ----------")
        else:
            res.append("---------- Alignment ----------")
        res.append("Size       : " + str(self.size))
        res.append("Type       : " + self.alignType)
        res.append("Score      : " + str(self.alignScore))
        if not self.isMultiple:
            res.append("Identity   : " + str(self.identity) \
                       + " ({0:.2f}%)".format(100 * (self.identity / len(self.seqA))))
            res.append("Similarity : " + str(self.similarity) \
                       + " ({0:.2f}%)".format(100 * (self.similarity / len(self.seqA))))
        res.append("Gaps       : " + (str(self.gaps) if self.gaps > 0 else "None"))
        res.append("")
        if self.isMultiple:
            res.append("PSSM : " + self.pssmDescription)
            res.append("Aligned seq. : " + self.seqA.getDescription())
            res.append("\t" + str(self.seqAGaps) + " Gaps, " + str(self.size - self.seqAGaps) +
                       " AAs (positions " + str(self.seqAStart) + " to " + str(self.seqAEnd) + ")")
        else:
            res.append("Upper seq. : " + self.seqA.getDescription())
            res.append("\t" + str(self.seqAGaps) + " Gaps, " + str(self.size - self.seqAGaps) +
                       " AAs (positions " + str(self.seqAStart) + " to " + str(self.seqAEnd) + ")")
            res.append("Lower seq. : " + self.seqB.getDescription())
            res.append("\t" + str(self.seqBGaps) + " Gaps, " + str(self.size - self.seqBGaps) +
                       " AAs (positions " + str(self.seqBStart) + " to " + str(self.seqBEnd) + ")")
        res.append("")

        if self.condensed:
            return "\n".join(res)

        listA, listB, listI = [], [], []
        maxALen, maxBLen = 0, 0
        seqAPos, seqBPos = self.seqAStart, self.seqBStart
        seqANextPos, seqBNextPos = seqAPos, seqBPos

        for index in range(self.size):
            a = self.seqA[index]
            if not a.isGap():
                seqANextPos += 1
            listA.append(str(a))
            maxALen = max(len(str(a)), maxALen)

            if not self.isMultiple:
                listI.append(self.seqInter[index])
                b = self.seqB[index]
                if not b.isGap():
                    seqBNextPos += 1
                listB.append(str(b))
                maxBLen = max(len(str(b)), maxBLen)

            if index != 0 and (index % self.chunkSize == 0 or index == len(self.seqA) - 1):
                res.append(str(seqAPos))
                for i in range(maxALen - 1, -1, -1):
                    res.append("".join([s[-i - 1] if len(s) > i else " " for s in listA]))

                if not self.isMultiple:
                    res.append("".join(listI))
                    for i in range(maxBLen):
                        res.append("".join([s[i] if len(s) > i else " " for s in listB]))
                        res.append(str(seqBPos))
                res.append("\n")
                listA, listB, listI = [], [], []
                maxALen, maxBLen = 0, 0
                seqAPos, seqBPos = seqANextPos, seqBNextPos

        return "\n".join(res)


class Align:
    """
    Represents an alignment matrix, used to determine the alignment score between two Sequences
    """

    def __init__(self, scoreMatrix):
        """
        Creates an AlignMatrix object that uses scoreMatrix as a scoring system between amino acids.
        """
        self._scoreMatrix = scoreMatrix  # scoring matrix to use

        self._colSeq = None  # copies of Sequences to align
        self._rowSeq = None

        self._alignMatrix = None  # Used for alignment scores between Sequences
        self._rowGapMatrix = None  # used for affine gap penalty
        self._colGapMatrix = None
        self._originMatrix = None  # Used for backtracking purposes

        self._alignMode = ""  # Align mode (global, semiglobal, local)
        self._isSuboptimal = False  # Is the alignment subobtimal ?
        self._isMultiple = False  # Is the alignment done against multiple Sequences at once ?
        self._resultCount = None  # The number of expected results
        self._subOptimalDepth = None

        self._maxAlignScore = 0  # Maximum score
        self._currentAlignScore = 0  # Score for this alignment

        self._maxScoreRows = []  # Index of maximum values from _alignMatrix
        self._maxScoreCols = []
        self._currentAlignPath = []
        self._bestAlignPath = []  # Sequence of origins for alignment with best score
        self._allAlignPaths = []  # all Sequences of origins for alignments with best scores

        self._iniGapPenalty = 0  # Initial gap penalty
        self._extGapPenalty = 0  # Extended gap penalty

        self._alignedRowSeq = None  # Aligned Sequences (result of alignment)
        self._alignedColSeq = None

    def globalAlign(self, seqA, seqB, iniGapPenalty=1, extGapPenalty=None, semiGlobal=False, resultCount=-1):
        """
        Returns all best global alignment between 'seqA' and 'seqB' with the provided
        initial and extended gap penalties.
        If 'semiGlobal' is True, allows alignment to discard the start and end of either Sequence.
        (this allows for better alignments when Sequences overlap only partially)
        """
        if semiGlobal:
            self._alignMode = "semiglobal"
        else:
            self._alignMode = "global"
        self._resultCount = resultCount

        # Initialize all data structures
        self.__initialize(seqA, seqB, iniGapPenalty, extGapPenalty)

        # Align in global mode, from the bottom right
        if self._alignMode == "global":
            self._currentAlignScore = self._alignMatrix[-1][-1]
            yield from self.__align(len(self._rowSeq), len(self._colSeq))
        else:
            self._currentAlignScore = self._maxAlignScore
            for row, col in zip(self._maxScoreRows, self._maxScoreCols):
                yield from self.__align(row, col)

    def localAlign(self, seqA, seqB, iniGapPenalty=1, extGapPenalty=None, resultCount=1, subOptimalDepth=0):
        """
        Returns n best local alignment between 'seqA' and 'seqB' with the provided
        initial and extended gap penalties, where n equals 'resultCount' (-1 for all).
        If 'subOptimalDepth' equals m, looks for m best suboptimal alignments as well.
        """
        self._alignMode = "local"

        # Initialize all data structures
        self.__initialize(seqA, seqB, iniGapPenalty, extGapPenalty)

        # Number of expected results for each max. score
        self._resultCount = resultCount

        # Align in local mode, from the maximum value (optimal)
        self._currentAlignScore = self._maxAlignScore
        for row, col in zip(self._maxScoreRows, self._maxScoreCols):
            yield from self.__align(row, col)  # Results

        # If suboptimal alignments are requested
        self._subOptimalDepth = 0
        while subOptimalDepth != 0:
            self._isSuboptimal = True
            subOptimalDepth -= 1
            self._subOptimalDepth += 1

            self.__clearBestPath()  # Reevalute scores with best alignment reset to 0

            # Number of expected results for each max. score
            self._resultCount = resultCount

            # Align in local mode, from the new maximum value (sub-optimal)
            self._currentAlignScore = self._maxAlignScore
            for row, col in zip(self._maxScoreRows, self._maxScoreCols):
                yield from self.__align(row, col)

    def multiAlign(self, Sequence, resultCount=1, subOptimalDepth=3):
        """
        Returns n best local M.S.A. of 'Sequence' against the ScoreMatrix (which must be 
        position-specific and provide gap penalties) where n equals 'resultCount' (-1 for all).
        If 'subOptimalDepth' equals m, looks for m best suboptimal alignments as well.
        """
        self._isMultiple = True
        colSeq = [i for i in range(len(self._scoreMatrix))]
        yield from self.localAlign(colSeq, Sequence, 0, None, resultCount, subOptimalDepth)

    def __initialize(self, seqA, seqB, iniGapPenalty, extGapPenalty):
        """
        Sets all initial values for required data structures.
        """
        # Sequences
        if len(seqA) == 0 or len(seqB) == 0:
            raise ValueError("Sequences to align cannot be empty")
        self._colSeq = seqA
        self._rowSeq = seqB

        # Alignments
        self._isSuboptimal = False
        self._subOptimalDepth = 0

        self._alignedRowSeq = Sequence()
        if self._isMultiple:
            # In MSA (align with PSSM) there is no column Sequence
            self._alignedColSeq = []
        else:
            self._alignedColSeq = Sequence()

        self._maxAlignScore = 0  # Maximum score found while aligning
        self._maxScoreRows = []  # Rows of maximum score
        self._maxScoreCols = []  # Columns of maximum score

        self._currentAlignPath = []  # indexes of the current alignment
        self._bestAlignPath = []  # indexes of the alignment with the best score

        # Gap penalties
        self._iniGapPenalty = iniGapPenalty
        self._extGapPenalty = iniGapPenalty if extGapPenalty is None else extGapPenalty

        # Matrices
        self._alignMatrix = [[0 for i in range(len(self._colSeq) + 1)] \
                             for j in range(len(self._rowSeq) + 1)]
        self._rowGapMatrix = deepcopy(self._alignMatrix)
        self._colGapMatrix = deepcopy(self._alignMatrix)

        self._originMatrix = [["" for i in range(len(self._colSeq) + 1)] \
                              for j in range(len(self._rowSeq) + 1)]

        # Global alignment : first line and colunm have initial scores and origins
        if self._alignMode == "global":
            self.__initAlignValues()

        # Fill all matrices
        for row in range(1, len(self._rowSeq) + 1):
            for col in range(1, len(self._colSeq) + 1):
                self.__fill(row, col)

        # Find best scores
        self.__findBestScore()

    def __initAlignValues(self):
        for i in range(1, max(len(self._colSeq), len(self._rowSeq)) + 1):
            val = self._iniGapPenalty + self._extGapPenalty * (i - 1)
            if i <= len(self._colSeq):  # First row
                self._alignMatrix[0][i] = val
                self._rowGapMatrix[0][i] = val
                self._colGapMatrix[0][i] = val
                self._originMatrix[0][i] = "L"  # Origin is left

            if i <= len(self._rowSeq):  # First column
                self._alignMatrix[i][0] = val
                self._rowGapMatrix[i][0] = val
                self._colGapMatrix[i][0] = val
                self._originMatrix[i][0] = "T"  # Origin is top

    def __fill(self, row, col):
        """
        Calculates the score for matrix values at row 'row' and column 'col'.
        """
        # Row gap matrix
        if self._isMultiple:
            self._rowGapMatrix[row][col] = self._scoreMatrix.getGapPenalty(col)
        else:
            initialGapScore = self._alignMatrix[row - 1][col] + self._iniGapPenalty
            extendedGapScore = self._rowGapMatrix[row - 1][col] + self._extGapPenalty
            self._rowGapMatrix[row][col] = max(initialGapScore, extendedGapScore)

        # Column gap matrix
        if self._isMultiple:
            self._colGapMatrix[row][col] = self._scoreMatrix.getGapPenalty(col - 1)
        else:
            initialGapScore = self._alignMatrix[row][col - 1] + self._iniGapPenalty
            extendedGapScore = self._colGapMatrix[row][col - 1] + self._extGapPenalty
            self._colGapMatrix[row][col] = max(initialGapScore, extendedGapScore)

        # Align matrix
        deleteScore = self._rowGapMatrix[row][col]
        insertScore = self._colGapMatrix[row][col]
        matchScore = self._alignMatrix[row - 1][col - 1] + \
                     self._scoreMatrix.getScore(self._rowSeq[row - 1], self._colSeq[col - 1])

        maxScore = 0
        if self._alignMode == "local":  # Local alignment
            self._alignMatrix[row][col] = maxScore = max(insertScore, deleteScore, matchScore, 0)
        else:  # Global alignment
            self._alignMatrix[row][col] = maxScore = max(insertScore, deleteScore, matchScore)

        # Origin matrix
        if insertScore == maxScore:
            self._originMatrix[row][col] += "L"
        if matchScore == maxScore:
            self._originMatrix[row][col] += "D"
        if deleteScore == maxScore:
            self._originMatrix[row][col] += "T"

    def __clearBestPath(self):
        """
        Resets all scores from the first best alignment and then reevaluates the affected scores.
        A new maximum score will be found from these, allowing for new alignment lookups.
        """
        self._allAlignPaths.extend(self._bestAlignPath)

        # Clear scores from best align path
        for row, col in self._bestAlignPath:
            self._alignMatrix[row][col] = 0  # Reset scores
            self._colGapMatrix[row][col] = 0
            self._rowGapMatrix[row][col] = 0
            self._originMatrix[row][col] = ""

        # Reevaluate scores for further rows and columns (affected values)
        row, col = self._bestAlignPath[-1]
        for furtherRow in range(row, len(self._rowSeq) + 1):
            for furtherCol in range(col, len(self._colSeq) + 1):
                if (furtherRow, furtherCol) not in self._allAlignPaths:
                    self._originMatrix[furtherRow][furtherCol] = ""
                    self.__fill(furtherRow, furtherCol)

        self._bestAlignPath = []
        self.__findBestScore()

    def __findBestScore(self):
        self._maxAlignScore = -1
        self._maxScoreCols = []
        self._maxScoreRows = []
        for row in range(1, len(self._rowSeq) + 1):
            for col in range(1, len(self._colSeq) + 1):
                score = self._alignMatrix[row][col]
                if score > self._maxAlignScore:  # New max alignment
                    self._maxAlignScore = score
                    self._maxScoreRows = [row]
                    self._maxScoreCols = [col]
                elif score == self._maxAlignScore:  # Same max alignment
                    self._maxScoreRows.append(row)
                    self._maxScoreCols.append(col)

    def __align(self, i, j):
        """
        Yields all best alignments starting from row i and column j.
        """
        if self._resultCount == 0:
            return

        # Local alignments are complete when we reach a null score
        # Global alignments are complete when we reach the beginning of the matrix
        # Semiglobal alignments are complete when we reach an edge of the matrix
        if (self._alignMode == "local" and self._alignMatrix[i][j] == 0) \
                or (self._alignMode == "global" and i == 0 and j == 0) \
                or (self._alignMode == "semiglobal" and (i == 0 or j == 0)):

            alignDescription = self._alignMode + \
                               ("-suboptimal" + "(" + str(self._subOptimalDepth) + ")") * self._isSuboptimal

            # Create result (Sequence obj. for MSA, Aligned obj. otherwise)
            if self._isMultiple:
                result = Aligned(self._alignedColSeq,
                                 Sequence(self._alignedRowSeq, self._rowSeq.getDescription()), j, i,
                                 alignDescription, self._currentAlignScore, self._scoreMatrix, True)
            else:
                result = Aligned(Sequence(self._alignedColSeq, self._colSeq.getDescription()),
                                 Sequence(self._alignedRowSeq, self._rowSeq.getDescription()), j, i,
                                 alignDescription, self._currentAlignScore, self._scoreMatrix, False)

            # Remember first best alignment for subobtimal lookup
            if self._bestAlignPath == []:
                self._bestAlignPath = deepcopy(self._currentAlignPath)

            self._resultCount -= 1
            yield result

        else:
            for origin in self._originMatrix[i][j]:
                self._currentAlignPath.append((i, j))

                if origin == "T":  # top
                    self._alignedColSeq.insert(0, "-")
                    self._alignedRowSeq.insert(0, self._rowSeq[i - 1])
                    yield from self.__align(i - 1, j)

                elif origin == "D":  # diagonal
                    self._alignedColSeq.insert(0, self._colSeq[j - 1])
                    self._alignedRowSeq.insert(0, self._rowSeq[i - 1])
                    yield from self.__align(i - 1, j - 1)

                elif origin == "L":  # left
                    self._alignedColSeq.insert(0, self._colSeq[j - 1])
                    self._alignedRowSeq.insert(0, "-")
                    yield from self.__align(i, j - 1)

                else:
                    raise ValueError("Origin must be T (top), D (diagonal) or L (left)")

                self._currentAlignPath.pop()
                del self._alignedColSeq[0]
                del self._alignedRowSeq[0]

    # Representation
    def __repr__(self):
        """
        Representation of the alignment matrix.
        """
        if self._colSeq is None:
            return "No alignment done yet"

        result = []
        sepSize = 5
        rep = {"": "", "L": "_   ", "D": " \\  ", "T": "   |", "LD": "_\\  ", "DT": " \\ |", "LT": "_  |",
               "LDT": "_\\ |"}
        # Amino Acid column Sequence
        result.append(" " * 2 * sepSize)
        for aa in self._colSeq:
            result[-1] += '{a!s:<{w}}'.format(a=aa, w=sepSize)

        # First line of values (no amino acid)
        result.append(" ")
        for origin in self._originMatrix[0]:
            result[-1] += '{o:>{w}}'.format(o=rep[origin], w=sepSize)

        result.append(" " * sepSize)
        for value in self._alignMatrix[0]:
            result[-1] += '{v:<{w}}'.format(v=value, w=sepSize)

        # Amino Acid line Sequence and values
        for values, origins, aa in zip(self._alignMatrix[1:], self._originMatrix[1:], self._rowSeq):
            result.append(" ")
            for origin in origins:
                result[-1] += '{o:>{w}}'.format(o=rep[origin], w=sepSize)
            result.append("")
            result[-1] += '{a!s:<{w}}'.format(a=aa, w=sepSize)
            for value in values:
                result[-1] += '{v:<{w}}'.format(v=value, w=sepSize)

        return "\n".join(result)
