# PyProt

**PyProt** (short for ***Python Proteins***) is a set of libraries designed to represent and manipulate proteins.

Here is an overview of the libraries and the tools they contain:
- `aminoacid` provides the `AminoAcid` class, used to represent any amino acid (or lack thereof)
- `sequence` provides the `Sequence` class, which represents an ordered sequence of amino acids (i.e. protein)
- `score` provides classes that attribute a score to substitutions between amino acids
- `align` provides classes that align two or more proteins together and represent these alignments in multiple ways
- `structure` implements algorithms that can be trained on datasets to issue structure predictions for new proteins
