from copy import deepcopy

from pyprot.base.aminoacid import AminoAcid


class Sequence(list):
    """
    Represents a sequence of amino acids.
    Inherits from list, and ensures all items are of type AminoAcid.
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

        # if copy constructor, copy attributes
        if isinstance(aminoAcids, Sequence):
            self._nameMode = aminoAcids._nameMode
            self._separator = aminoAcids._separator
            if self._description == "":
                self._description = aminoAcids._description

        # format aminoAcids into a list of AminoAcid objects, and add it to list
        self.extend(Sequence.__formatList(aminoAcids))

    @staticmethod
    def __formatList(aminoAcids):
        """Formats 'aminoAcids' into a list of AminoAcid objects."""

        # No amino acids result in an empty list
        if aminoAcids is None:
            return []

        # A single AminoAcid is copied and put within a list
        elif isinstance(aminoAcids, AminoAcid):
            return [AminoAcid(aminoAcids)]  # Copy constructor

        # A string is converted to a list, based on its
        elif isinstance(aminoAcids, str):
            if aminoAcids.isupper():  # Multiple Amino Acids in short name mode
                return [AminoAcid(aa) for aa in aminoAcids]
            else:  # A single Amino Acid with any name mode
                return [AminoAcid(aminoAcids)]

        # A list is copied with all of its items converted to AminoAcid objects
        elif isinstance(aminoAcids, list):
            return [AminoAcid(aa) for aa in aminoAcids]

        # No other supported types
        else:
            raise TypeError("aminoAcids must be a Sequence, list, AminoAcid object, string or None")

    def __repr__(self):
        """Representation"""
        return str(self)

    def __str__(self):
        """String conversion"""
        return self._separator.join([aa.getName(self._nameMode) for aa in self])

    def setDescription(self, description):
        """Sets the base's description"""
        self._description = description

    def getDescription(self):
        """Returns the base's description."""
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

    def __setitem__(self, key, value):
        """Sets value for a slice of the sequence"""
        list.__setitem__(self, key, AminoAcid(value))

    def insert(self, index, aminoAcids):
        """
        Inserts aminoAcids into the base at index 'index'.
        List objects will not be embedded as is, instead their items will be inserted in the same order, individually.
        @param aminoAcids must be compatible with the Sequence constructor
        @param index is the index at which aminoAcids is inserted
        """
        for aa in Sequence.__formatList(aminoAcids):
            list.insert(self, index, aa)
            index += 1

    def extend(self, aminoAcids):
        """
        Extends the base by adding 'aminoAcids' at its end.
        @param aminoAcids must be compatible with the Sequence constructor
        """
        list.extend(self, Sequence.__formatList(aminoAcids))
