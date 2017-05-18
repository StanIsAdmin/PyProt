from pyprot.base.protein import Protein


def getProteinsFromFasta(path):
    """
    Loads the FASTA file located in 'path' and yields the Proteins it contains.
    """
    with open(path, 'r') as fastaFile:
        newProtein = None
        for line in fastaFile:
            line_s = line.strip()
            if line_s != "" and line_s[0] == ">":
                if newProtein is not None:
                    yield newProtein
                newProtein = Protein(None, line_s[1:])
            else:
                newProtein.extend(line_s)
        if len(newProtein) > 0:
            yield newProtein
