# PyProt

## What's in the box
**PyProt** (Short for _Python Proteins_) is a python package with libraries designed to represent and manipulate proteins.

Here is an overview of the sub-packages and the libraries they contain:
- `base` contains basic representation classes
  - `aminoacid` defines the `AminoAcid` class which represents a single amino acid
  - `sequence` represents an amino acid `Sequence` (or _protein_) and works as a python `list`
- `data` contains parsers for standard data files
  - `fasta` parses and saves `.fasta` files which contain proteins
  - `dssp` parses and saves `.dssp` files which contain score matrices
- `align` contains classes that align proteins together and represent the results
  - `align` defines the `Align` class which aligns two or more sequences together
  - `aligned` stores the results of the alignment in the `Aligned` class
  - `blosum` creates scoring matrices with the `BLOSUM` algorithm
  - `score` represents scoring matrices, both position-specific (`PSSM`) and not (`ScoreMatrix`)
- `structure` implements algorithms that can be trained on data sets to issue structure predictions for new proteins
  - `GOR` implements the `GOR.3` algorithm for structure prediction

## Dependencies
Python 3

## Installation
You can install this package in two different ways:

With `pip`, by running

`pip install git+https://github.com/jkbr/httpie.git`

By cloning this repo and adding the repo to your python `PATH`


## Testing
Unit tests are embedded in the corresponding packages.

## Contribute
Anyone is welcome to contribute by submitting a [pull request](https://help.github.com/articles/about-pull-requests/) or opening [new issues](https://help.github.com/articles/about-issues/).