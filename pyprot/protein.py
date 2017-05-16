from copy import deepcopy

aminoAcidNames = (
    ("none/gap", "gap", "-"),
    ("alanine", "ala", "A"),
    ("cysteine", "cys", "C"),
    ("aspartic acid", "asp", "D"),
    ("glutamic acid", "glu", "E"),
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
    ("asparagine/aspartic acid", "asx", "B"),
    ("glutamine/glutamic acid", "glx", "Z"),
    ("leucine/isoleucine", "xle", "J"),
    ("selenocysteine", "sec", "U"),
    ("pyrrolysine", "pyl", "O"),
    ("undetermined", "xaa", "X")
)


class AminoAcid:
    """
    Represents one of the amino acids that can be found in genetic sequences.
    Can be one of the twenty-two amino acids, four undetermined combinations of possible amino acids, and gaps.
    """

    # Dictionary mapping name to id
    _nameDict = {aminoAcidNames[id][i]: id for i in range(3) for id in range(len(aminoAcidNames))}

    _nameModes = {"long": 0, "medium": 1, "short": 2}  # choices for name length
    _defaultNameMode = "short"  # short name by default

    def __init__(self, aminoAcid):
        """
        Creates an AminoAcid object representing one of the possible Amino Acids.
        aminoAcid can be the name of an amino acid, or an AminoAcid object (in which case a copy is created).
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
    def getAllNames(nameMode=_defaultNameMode):
        try:
            nameIndex = AminoAcid._nameModes[nameMode]  # get index of name mode
        except:
            raise TypeError("nameMode must be 'short', 'medium' or 'long'")

        for aa in aminoAcidNames[1:]:  # we exclude the gap (first item)
            yield aa[nameIndex]

    # Representation
    def __repr__(self):
        nameIndex = AminoAcid._nameModes[self._defaultNameMode]
        return aminoAcidNames[self._id][nameIndex]  # default name mode

    def __str__(self):
        return self.getName()  # default name mode

    def getName(self, nameMode=_defaultNameMode):
        try:
            nameIndex = AminoAcid._nameModes[nameMode]  # get index of name mode
        except:
            raise TypeError("nameMode must be 'short', 'medium' or 'long'")

        return aminoAcidNames[self._id][nameIndex]

    # Comparison
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

    def isGap(self):
        return self._id == 0

    # Hashing
    def __hash__(self):
        return hash(self._id)


class Sequence:
    """
    Represents a sequence of amino acids.
    """

    def __init__(self, aminoAcids=None, description=""):
        """
        Creates a Sequence object that represents the amino acid sequence contained in aminoAcids.
        aminoAcids can be one of the following :
        - None, meaning the Sequence is empty (default)
        - an AminoAcid object
        - a string of X AminoAcid short (uppercase) names or 1 AminoAcid name
        - a list containing AminoAcid objects and/or strings of individual AminoAcid names
        """

        self._nameMode = "short"  # the way in which AA names are displayed
        self._separator = ""  # how to separate AA names when displayed
        self._description = description  # description of the sequence
        if self._description == "" and isinstance(aminoAcids, Sequence):
            self._description = aminoAcids.getDescription()

        # Format aminoAcids into a list of AminoAcid objects.
        self._aaList = self.__formatAAList(aminoAcids)  # List of amino acids

    @staticmethod
    def __formatAAList(aminoAcids):
        """Formats 'aminoAcids' into a list of AminoAcid objects."""

        # Sequence object's aaList is deep copied
        if isinstance(aminoAcids, Sequence):
            return deepcopy(aminoAcids._aaList)

        # None becomes an empty Sequence
        elif aminoAcids is None:
            return []

        # A string is converted to a list
        elif isinstance(aminoAcids, str):
            if aminoAcids.isupper():  # Multiple Amino Acids in short name mode
                return [AminoAcid(aa) for aa in aminoAcids]
            else:  # A single Amino Acid with any name mode
                return [AminoAcid(aminoAcids)]

        # An AminoAcid is copied (from name) and put within a list
        elif isinstance(aminoAcids, AminoAcid):
            return [AminoAcid(aminoAcids)]  # Copy constructor

        # A list is copied with all of its items converted to AminoAcid objects
        elif isinstance(aminoAcids, list):
            return [AminoAcid(aa) for aa in aminoAcids]

        # No more supported types
        else:
            raise TypeError("aminoAcids must be a Sequence, list, AminoAcid object, string or None")

    # Size and comparison
    def __len__(self):
        return len(self._aaList)

    def __gt__(self, other):
        return len(self) > len(other)

    def __lt__(self, other):
        return len(self) < len(other)

    def __ge__(self, other):
        return len(self) >= len(other)

    def __le__(self, other):
        return len(self) <= len(other)

    def __eq__(self, other):
        return len(self) == len(other) and all(a == b for a, b in zip(self, other))

    def __ne__(self, other):
        return not self == other

    # Iteration
    def __iter__(self):
        return iter(self._aaList)

    # Representation
    def __repr__(self):
        """Representation"""
        return str(self)

    def __str__(self):
        """String conversion"""
        return self._separator.join([aa.getName(self._nameMode) for aa in self])

    def getDescription(self):
        """Returns the description of the sequence."""
        return self._description

    def setNameMode(self, newMode):
        """Changes the name display mode to 'newMode'."""
        if newMode in ("long", "medium", "short"):
            self._nameMode = newMode
        else:
            raise ValueError("newMode must be 'long', 'medium' or 'short'")

    def setSeparator(self, newSep):
        """Changes the string that separates each displayed AminoAcid."""
        self._separator = newSep

    # Item and slice manipulation
    def __getitem__(self, key):
        """Return a Sequence object containing copies of the items from the slice"""
        if isinstance(key, slice):
            return Sequence([self._aaList[index] for index in range(*key.indices(len(self)))])
        else:
            try:
                return self._aaList[key]  # Return the item of index 'key'
            except:
                raise ValueError("key does not represent an index or slice")

    def __setitem__(self, key, value):
        """Set value for a slice of the sequence"""
        if isinstance(key, slice):
            for index in range(*key.indices(len(self))):  # range(start, stop, step)
                self._aaList[index] = AminoAcid(value)  # create copies
        else:
            try:
                self._aaList[key] = AminoAcid(value)  # Set value for one item
            except:
                raise ValueError("key does not represent an index or slice")

    def __delitem__(self, key):
        """Delete a slice of the sequence"""
        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            del self._aaList[start:stop:step]
        else:
            try:
                del self._aaList[key]  # Delete one item
            except:
                raise ValueError("key does not represent an index or slice")

    # Sequence modification
    def insert(self, index, subSequence):
        """
        Inserts subSequence into the Sequence at index 'index' (default is 0).
        subSequence must be compatible with an AminoAcid list, as specified by __formatAAList.
        """
        # We need a formatted sequence
        for aa in self.__formatAAList(subSequence):
            self._aaList.insert(index, aa)
            index += 1

    def extend(self, subSequence):
        """
        Same as calling insert at the end of the Sequence.
        """
        self.insert(len(self), subSequence)

    def remove(self, subSequence, startIndex=0):
        """
        Removes the first occurrence of 'subSequence' in self, starting at index. Returns
        - start index of the sub-sequence within self (if it exists)
        - False if subSequence is not a sub-sequence of self
        """
        lookupResult = self.hasSubSequence(subSequence, startIndex)
        if isinstance(lookupResult, int):
            del self[startIndex:startIndex + len(subSequence)]

        return lookupResult

    def delete(self, start=0, stop=None, step=None):
        """
        Deletes Amino Acids between indexes start (included) and stop (excluded).
        If start is not specified, deletion will begin at index 0.
        If stop is not specified, deletion will stop after one item.
        """
        if stop is None:
            stop = start + 1
        if step is None:
            step = 1
        del self[start:stop:step]

    # Lookup
    def __contains__(self, item):
        """
        Returns True if item is contained in Sequence, False if not.
        """
        return item in self._aaList

    def hasSubSequence(self, subSequence, startIndex=0):
        """
        Checks if 'subSequence' is a sub-sequence of self, starting at startIndex. Returns
        - start index of the sub-sequence within self, if it exists
        - False if sequence is not a sub-sequence of self
        """
        subLen = len(subSequence)
        endIndex = len(self) - subLen + 1
        if endIndex < startIndex or startIndex < 0:
            raise ValueError("subSequence does not fit in sequence after startIndex")

        for i in range(startIndex, endIndex):
            if all(self[i + j] == subSequence[j] for j in range(subLen)):
                return i

        return False


def loadFasta(path):
    """
    Loads the FASTA file located in 'path' and yields the Sequences it contains.
    """
    with open(path, 'r') as fastaFile:
        newSequence = None
        for line in fastaFile:
            line_s = line.strip()
            if line_s != "" and line_s[0] == ">":
                if newSequence is not None:
                    yield newSequence
                newSequence = Sequence(None, line_s[1:])
            else:
                newSequence.extend(line_s)
        if len(newSequence) > 0:
            yield newSequence
