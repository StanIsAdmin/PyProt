from copy import deepcopy

AA_NAMES = (
    ("alanine", "ala", "A"),
    ("cysteine", "cys", "C"),
    ("aspartate", "asp", "D"),
    ("glutamate", "glu", "E"),
    ("phenylalanine", "phe", "F"),
    ("glycine", "gly", "G"),
    ("histidine", "his", "H"),
    ("isoleucine", "ile", "I"),
    ("lysine", "lys", "K"),
    ("leucine", "leu", "L"),
    ("methionine", "met", "M"),
    ("asparagine", "asn", "N"),
    ("proline", "pro", "P"),
    ("glutamine", "gln", "Q"),
    ("arginine", "arg", "R"),
    ("serine", "ser", "S"),
    ("threonine", "thr", "T"),
    ("valine", "val", "V"),
    ("tryptophan", "trp", "W"),
    ("tyrosine", "tyr", "Y"),
    ("selenocysteine", "sec", "U"),
    ("pyrrolysine", "pyl", "O"),
    ("asparagine/aspartate", "asx", "B"),
    ("glutamine/glutamate", "glx", "Z"),
    ("leucine/isoleucine", "xle", "J"),
    ("undetermined", "xaa", "X"),
    ("gap", "gap", "-"),
    ("termination", "term", "|")
)

AA_NAMES_CLASSIC_RANGE = (0, 20)
AA_NAMES_EXTENDED_RANGE = (0, 26)
AA_NAMES_GAP_INDEX = 26
AA_NAMES_TERM_INDEX = 27


class AminoAcid:
    """
    Represents one of the amino acids that can be found in genetic sequences.
    Can be one of the following :
    - any of the twenty amino acids
    - any of four combinations of possible amino acids
    - selenocysteine, pyrrolysine, a gap or termination codon
    The full list of possible amino acids is defined by AA_NAMES.
    """

    # Dictionary mapping name to id
    _nameDict = {AA_NAMES[id][i]: id for i in range(3) for id in range(len(AA_NAMES))}
    
    _nameModes = {"long": 0, "medium": 1, "short": 2}  # choices for name length
    _defaultNameMode = "short"  # short name by default

    def __init__(self, aminoAcid):
        """
        Creates an AminoAcid object representing one of the possible amino acids.
        @param aminoAcid can be the name of an amino acid, or an AminoAcid object (in which case a copy is created).
        """
        self._id = None  # id of the amino acid within the name group

        if isinstance(aminoAcid, str):
            if len(aminoAcid) == 1:
                self._id = self.__getIdByName(aminoAcid.upper())  # id from short aminoAcid name
            else:
                self._id = self.__getIdByName(aminoAcid.lower())  # id from other aminoAcid name
        elif isinstance(aminoAcid, AminoAcid):
            self._id = aminoAcid._id  # copy of id
        else:
            raise TypeError("aminoAcid must be a string or an AminoAcid object")

    @staticmethod
    def __getIdByName(name):
        try:
            return AminoAcid._nameDict[name]  # get index of name mode
        except:
            raise ValueError("Could not find amino acid name {}".format(name))

    @staticmethod
    def __getNameModeIndex(nameMode):
        try:
            return AminoAcid._nameModes[nameMode]  # get index of name mode
        except:
            raise TypeError("nameMode must be 'short', 'medium' or 'long'")

    @staticmethod
    def getNames(nameMode=_defaultNameMode):
        """Yields the names of the 20 basic amino acids."""
        start, stop = AA_NAMES_CLASSIC_RANGE
        yield from AminoAcid.getNamesInRange(start, stop, nameMode)

    @staticmethod
    def getAllNames(nameMode=_defaultNameMode):
        """Yields the names of all represented amino acids, excepting gaps and termination codons."""
        start, stop = AA_NAMES_EXTENDED_RANGE
        yield from AminoAcid.getNamesInRange(start, stop, nameMode)

    @staticmethod
    def getNamesInRange(startIndex, stopIndex, nameMode=_defaultNameMode):
        """Yields the names of amino acids in AA_NAMES, from startIndex to stopIndex (excluded)."""
        nameModeIndex = AminoAcid.__getNameModeIndex(nameMode)
        for aa in AA_NAMES[startIndex:stopIndex]:
            yield aa[nameModeIndex]

    def isGap(self):
        """True if this amino acid is a gap, false otherwise."""
        return self._id == AA_NAMES_GAP_INDEX

    def isTermination(self):
        """True if this amino acid is a termination codon, false otherwise."""
        return self._id == AA_NAMES_TERM_INDEX

    def getName(self, nameMode=_defaultNameMode):
        try:
            nameIndex = AminoAcid._nameModes[nameMode]  # get index of name mode
        except:
            raise TypeError("nameMode must be 'short', 'medium' or 'long'")

        return AA_NAMES[self._id][nameIndex]

    def __repr__(self):
        """Equivalent to getName()"""
        return self.getName()

    def __str__(self):
        """Equivalent to getName()."""
        return self.getName()  # default name mode

    # Comparison and hashing allow to manipulate and sort instances more efficiently
    # these functions do not have any biological meaning and their results may change over time.
    def __eq__(self, other):
        return self._id == other._id

    def __ne__(self, other):
        return self._id != other._id

    def __gt__(self, other):
        return self._id > other._id

    def __ge__(self, other):
        return self._id >= other._id

    def __lt__(self, other):
        return self._id < other._id

    def __le__(self, other):
        return self._id <= other._id

    def __hash__(self):
        return hash(self._id)
