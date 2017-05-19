from pyprot.base.sequence import Sequence


def getSequencesFromFasta(path):
    """
    Loads the FASTA file located in 'path' and yields the Sequences it contains.
    """
    with open(path, 'r') as fastaFile:
        newProtein = None
        for line in fastaFile:
            line_s = line.strip()
            if line_s != "" and line_s[0] == ">":
                if newProtein is not None:
                    yield newProtein
                newProtein = Sequence(None, line_s[1:])
            else:
                newProtein.extend(line_s)
        if len(newProtein) > 0:
            yield newProtein
