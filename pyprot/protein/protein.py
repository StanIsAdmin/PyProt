from copy import deepcopy

from pyprot.protein.aminoacid import AminoAcid


class Protein:
    """
    Represents a sequence of amino acids.
    """

    def __init__(self, aminoAcids=None, description=""):
        """
        Creates a Protein object that represents the amino acid sequence contained in aminoAcids.
        aminoAcids can be one of the following :
        - None, meaning the Protein is empty (default)
        - an AminoAcid object
        - a string of X AminoAcid short (uppercase) names or 1 AminoAcid name
        - a list containing AminoAcid objects and/or strings of individual AminoAcid names
        """

        self._nameMode = "short"  # the way in which AA names are displayed
        self._separator = ""  # how to separate AA names when displayed
        self._description = description  # description of the sequence
        if self._description == "" and isinstance(aminoAcids, Protein):
            self._description = aminoAcids.getDescription()

        # Format aminoAcids into a list of AminoAcid objects.
        self._aaList = self.__formatAAList(aminoAcids)  # List of amino acids

    @staticmethod
    def __formatAAList(aminoAcids):
        """Formats 'aminoAcids' into a list of AminoAcid objects."""

        # Protein object's aaList is deep copied
        if isinstance(aminoAcids, Protein):
            return deepcopy(aminoAcids._aaList)

        # None becomes an empty Protein
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
            raise TypeError("aminoAcids must be a Protein, list, AminoAcid object, string or None")

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
        """Return a Protein object containing copies of the items from the slice"""
        if isinstance(key, slice):
            return Protein([self._aaList[index] for index in range(*key.indices(len(self)))])
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

    # Protein modification
    def insert(self, index, subSequence):
        """
        Inserts subSequence into the Protein at index 'index' (default is 0).
        subSequence must be compatible with an AminoAcid list, as specified by __formatAAList.
        """
        # We need a formatted sequence
        for aa in self.__formatAAList(subSequence):
            self._aaList.insert(index, aa)
            index += 1

    def extend(self, subSequence):
        """
        Same as calling insert at the end of the Protein.
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
        Returns True if item is contained in Protein, False if not.
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
