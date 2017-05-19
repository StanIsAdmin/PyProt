# PyProt

**PyProt** (short for _Python Proteins_) is a python package with libraries designed to represent and manipulate proteins.

## What's in the box

Here is an overview of the sub-packages and the libraries they contain:
- `base` contains basic representation classes
  - `aminoacid` defines the `AminoAcid` class which represents a single amino acid
  - `sequence` represents an amino acid `Sequence` (or _protein_) and works as a python `list`
- `data` contains parsers for standard data files
  - `fasta` parses and saves `.fasta` files which contain proteins
  - `dssp` parses and saves `.dssp` files which contain score matrices
- `align` contains classes that align proteins together and represent the results
  - `align` defines the `Align` class which aligns sequences together, and `Aligned` which stores the alignment results
  - `blosum` creates scoring matrices with the `BLOSUM` algorithm
  - `score` represents scoring matrices, both position-specific (`PSSM`) and not (`ScoreMatrix`)
- `structure` implements algorithms that can be trained on data sets to issue structure predictions for new proteins
  - `GOR` implements the `GOR.3` algorithm for structure prediction

## Dependencies
In order to use this package, you'll need a working version of [Python 3](https://www.python.org/download/releases/3.0/).

Installing the package with the following instructions will automatically install its dependencies :
- [matplotlib](https://matplotlib.org/)

## Installation
You can install this package in the following ways:

- By executing the following in your command line

`pip install git+https://github.com/StanIsAdmin/PyProt.git#egg=pyprot`

- By downloading the package's source code [here](https://github.com/StanIsAdmin/PyProt/archive/master.zip), unzipping it and then running

`pip install <downloaded-code-path>`

- Once the source code is downloaded and extracted, you can also copy-paste it to any directory in your Python `PATH`.

## Contribute
Anyone is welcome to contribute by submitting a [pull request](https://help.github.com/articles/about-pull-requests/) or by opening [new issues](https://help.github.com/articles/about-issues/).
