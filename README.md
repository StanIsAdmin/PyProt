# PyProt

**PyProt** (short for ***Python Proteins***) is a set of libraries designed to represent and manipulate proteins.

Here is an overview of the libraries and the tools they contain (classes in `CamelCase` and functions in `camelCase`):
- `protein` provides tools to represent amino acids as well as sequences of amino acids (i.e. proteins)
 - `AminoAcid` represents a single, possibly undetermined amino acid or gap from the `aminoAcidNames` tuple
 - `Sequence` represents a protein and can be manipulated as a python list
 - `loadFasta` parses `.fasta` files and generates the resulting `Sequence` objects
- `score` implements scoring systems for substitutions between amino acids
 - `ScoreMatrix` is a standard substitution-based score matrix
 - `PSSM` is a position-specific score matrix
 - `blosumFromFasta` creates a `ScoreMatrix` from a `.fasta` file based on the BLOSUM method
- `align` contains classes that align two or more proteins together and represent the results
 - `Aligned` is the result of any alignment, and contains information about the aligned sequences
 - `Align` is a class with multiple functions to satisfy all your alignment needs
- `structure` implements algorithms that can be trained on datasets to issue structure predictions for new proteins
 - `DSSP` can be used to parse `.dssp` files
 - `GOR3` implements the eponymous algorithm
