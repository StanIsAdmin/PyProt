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
                    res.append(str(seqBPos))
                    for i in range(maxBLen):
                        res.append("".join([s[i] if len(s) > i else " " for s in listB]))
                res.append("\n")
                listA, listB, listI = [], [], []
                maxALen, maxBLen = 0, 0
                seqAPos, seqBPos = seqANextPos, seqBNextPos

        return "\n".join(res)
