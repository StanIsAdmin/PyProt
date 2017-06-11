# PyProt

[![Build Status](https://travis-ci.org/StanIsAdmin/PyProt.svg?branch=master)](https://travis-ci.org/StanIsAdmin/PyProt)
[![Coverage Status](https://coveralls.io/repos/github/StanIsAdmin/PyProt/badge.svg?branch=master)](https://coveralls.io/github/StanIsAdmin/PyProt?branch=master)
[![PyPi Package](https://img.shields.io/badge/pypi%20package-0.1--alpha-orange.svg)](https://pypi.python.org/pypi/pyprot)



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
In order to use this package, you'll need a working version of [Python 3.3](https://www.python.org/download/releases/3.0/) or later installed, as well as [pip](https://pypi.python.org/pypi/pip).

The installation process will automatically install all of the package's dependencies, which are listed in the `setup-req.txt` file.

## Installation
You can install `Pyprot` in the following ways (make sure you use a `Python 3` version of `pip`):

- By executing the following in your command line

`pip install git+https://github.com/StanIsAdmin/PyProt.git --user`

- By downloading the package's source code [here](https://github.com/StanIsAdmin/PyProt/archive/master.zip), unzipping it and then running

`pip install <downloaded-code-path> --user`

## Documentation

Not sure what's in the box yet ? Check the [online documentation](https://stanisadmin.github.io/PyProt/).

The source code is documented in the standard `docstring` format, so its documentation will appear automatically if you use an editor that supports that format (which really means any editor but vim).

Working examples are provided in the `examples/` folder of the repository, with code and explanations embedded inside `jupyter notebook` files. You can read them from GitHub, but in order to run them yourself, you'll need to install [Jupyter](http://jupyter.org/).

## Contribute
Anyone is welcome to contribute by submitting a [pull request](https://help.github.com/articles/about-pull-requests/) or by opening [new issues](https://help.github.com/articles/about-issues/).
