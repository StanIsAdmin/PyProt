from math import log, sqrt

import matplotlib.pyplot as pyplot
import numpy as np

from pyprot.base.aminoacid import AminoAcid


class GOR3:
    """
    Implements the GOR III secondary structure prediction algorithm.
    Objects of this class must be trained with known examples of sequences and their structure
    before being able to predict the structures of new sequences.
    """

    def __init__(self):
        self.structures = "HETC"
        self.roc = {(s, r): 0 for s in self.structures for r in ("TP", "TN", "FP", "FN")}
        self.minScore, self.maxScore = 0, 0
        self.scores = {(classifier, realS): [] for classifier in self.structures for realS in self.structures}
        self.correctPred = 0
        self.totalPred = 0
        self.neighbourOffset = 8

        self.trainings = 0  # Number of trainings (one per AA)
        self.strucCount = {s: 0 for s in self.structures}
        self.pairCount = {(s, a): 0 for s in self.structures for a in AminoAcid.getAllNames()}
        self.tripletCount = {}
        for s in self.structures:
            for a in AminoAcid.getAllNames():
                for na in AminoAcid.getAllNames():  # Neighbour AA
                    self.tripletCount[(s, a, na)] = 0

    def train(self, sequence, structure):
        """
        Trains the system with a known example of a sequence and its structure.
        """
        self.trainings += len(sequence)

        for index in range(len(sequence)):
            curAminoacid = sequence[index]
            curStructure = structure[index]
            self.strucCount[curStructure] += 1
            self.pairCount[(curStructure, curAminoacid)] += 1

            for neiAminoacid in self.neighbourValues(sequence, index):
                self.tripletCount[(curStructure, curAminoacid, neiAminoacid)] += 1

    def predict(self, sequence, realStructure=None):
        """
        Returns the predicted structure of 'sequence', based on received training.
        """
        structure = []  # Result: predicted structure

        # Predict structures for each aminoacid in sequence
        for index in range(len(sequence)):
            curAminoacid = sequence[index]

            # First possible structure
            predStructure = self.structures[0]
            predScore = self.__getScore(sequence, index, predStructure)

            # Other structures
            for curStructure in self.structures[1:]:
                curScore = self.__getScore(sequence, index, curStructure)

                # Remember structure that gives best score
                if curScore > predScore:
                    predStructure = curStructure
                    predScore = curScore

            structure.append(predStructure)
        structure = "".join(structure)

        if not realStructure is None:
            self.__quality(sequence, structure, realStructure)

        return structure

    def __getScore(self, sequence, index, struct):
        """
        Returns I(deltaS, R) as defined by the GOR III algorithm.
        """
        aminoacid = sequence[index]
        scoreTerms = []

        score = self.pairCount[(struct, aminoacid)]
        score /= sum(self.pairCount[(otherStruct, aminoacid)] for otherStruct in self.getStructures(struct))
        scoreTerms.append(score)

        score = (self.trainings - self.strucCount[struct]) / self.strucCount[struct]
        scoreTerms.append(score)

        for neiAminoacid in self.neighbourValues(sequence, index):
            score = self.tripletCount[(struct, aminoacid, neiAminoacid)]
            score /= sum(
                self.tripletCount[(otherStruct, aminoacid, neiAminoacid)] for otherStruct in self.getStructures(struct))
            scoreTerms.append(score)

            score = sum(self.pairCount[(otherStruct, aminoacid)] for otherStruct in self.getStructures(struct))
            score /= self.pairCount[(struct, aminoacid)]
            scoreTerms.append(score)

        return sum(log(s) for s in scoreTerms)

    def neighbourOffsets(self):
        for offset in range(-self.neighbourOffset, self.neighbourOffset + 1):
            if offset != 0:
                yield offset

    def neighbourValues(self, sequence, index):
        for offset in self.neighbourOffsets():
            neiIndex = index + offset
            if neiIndex >= 0 and neiIndex < len(sequence):
                yield sequence[neiIndex]

    def getStructures(self, exclude=None):
        for s in self.structures:
            if s != exclude:
                yield s

    def __quality(self, sequence, structure, reality):
        self.totalPred += len(structure)  # Total predictions
        self.correctPred += sum([1 if s == r else 0 for s, r in zip(structure, reality)])  # Correct predictions
        for index in range(len(sequence)):
            scores = [self.__getScore(sequence, index, classifier) for classifier in self.structures]
            maxScore = max(scores)
            for classifier, score in zip(self.structures, scores):
                realS = reality[index]
                self.scores[(classifier, realS)].append(score)
                self.minScore = score if score < self.minScore else self.minScore
                self.maxScore = score if score > self.maxScore else self.maxScore
                if score == maxScore:  # Positive
                    if classifier == realS:  # True
                        self.roc[(classifier, "TP")] += 1
                    else:  # False
                        self.roc[(classifier, "FP")] += 1
                else:  # Negative
                    if classifier != realS:  # True
                        self.roc[(classifier, "TN")] += 1
                    else:  # False
                        self.roc[(classifier, "FN")] += 1

    def getQuality(self):
        tp = sum([self.roc[(c, "TP")] for c in self.structures])
        tn = sum([self.roc[(c, "TN")] for c in self.structures])
        fp = sum([self.roc[(c, "FP")] for c in self.structures])
        fn = sum([self.roc[(c, "FN")] for c in self.structures])
        mcc = (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        q3 = self.correctPred / self.totalPred
        return q3, mcc

    def plotROC(self):
        for classifier in self.structures:
            x, y = [], []
            for threshold in np.arange(self.minScore, self.maxScore, 0.1):
                tp, tn, fp, fn = 0, 0, 0, 0
                for realS in self.structures:
                    for score in self.scores[(classifier, realS)]:
                        if score >= threshold:  # Positive
                            if classifier == realS:  # True
                                tp += 1
                            else:  # False
                                fp += 1
                        else:  # Negative
                            if classifier != realS:  # True
                                tn += 1
                            else:  # False
                                fn += 1
                tpr = tp / (tp + fn)
                fpr = fp / (tn + fp)
                x.append(fpr)
                y.append(tpr)
            print("----- Receiver Operating Characteristic Curve -----")
            print("Classifier:", classifier)
            x.sort()
            y.sort()
            auc = np.trapz(y, x)  # Area Under Curve
            print("AUC:", auc)
            pyplot.title("Classifier " + classifier)
            pyplot.xlabel('False Positive Rate')
            pyplot.ylabel('True Positive Rate')
            pyplot.plot([0, 1], [0, 1])
            pyplot.plot(x, y)
            pyplot.show()
