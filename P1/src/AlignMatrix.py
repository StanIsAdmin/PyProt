class AlignMatrix:
    """
    Represents an alignment matrix, used to determine the alignment score between two sequences
    """
    
    def __init__(self, scoreMatrix, seqA, seqB):
        if not isinstance(seqA, Sequence) ans isinstance(seqB, Sequence):
            raise TypeError("seqA and seqB must be Sequence objects")
        if len(seqA)==0 or len(seqB)==0:
            raise ValueError("seqA and seqB cannot be empty")
            
        if not isinstance(scoreMatrix, Score):
            raise TypeError("scoreMatrix must be a Score object")
            
        self._seqA = seqA
        self._seqB = seqB
        self._scoreMatrix = scoreMatrix
        
        
    def 
        