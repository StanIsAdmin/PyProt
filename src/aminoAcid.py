from copy import deepcopy as deepcopy

nameGroups = (
		("none/gap", "gap", "-"),
		("alanine", "ala", "A"),
		("cysteine", "cys", "C"),
		("aspartic acid", "asp", "D"),
		("glutamic acid", "glu", "E"),
		("phenylalanine", "phe", "F"),
		("glycine", "gly", "G"),
		("histidine", "his", "H"),
		("isoleucine", "ile", "I"),
		("lysine", "lys", "K"),
		("leucine", "leu", "L"),
		("methionine", "met", "M"),
		("asparagine", "asn", "N"),
		("proline", "pro", "P"),
		("glutamine", "gln", "Q"),
		("arginine", "arg", "R"),
		("serine", "ser", "S"),
		("threonine", "thr", "T"),
		("valine", "val", "V"),
		("tryptophan", "trp", "W"),
		("tyrosine", "tyr", "Y"),
		("asparagine/aspartic acid", "asx", "B"),
		("glutamine/glutamic acid", "glx", "Z"),
		("leucine/isoleucine", "xle", "J"),
		("selenocysteine", "sec", "U"),
		("pyrrolysine", "pyl", "O"),
		("undetermined", "und", "X")
	)

class AminoAcid:
	"""
	Represents one of the amino acids that can be found in genetic sequences.
	Can be one of the twenty-two amino acids, four undetermined combinations of possible amino acids, and gaps.
	"""
	#Tuple of tuples listing amino acid names
	_nameGroups = deepcopy(nameGroups)
	
	#Dictionary mapping name to id
	_nameDict = {nameGroups[id][i]:id for i in range(3) for id in range(len(nameGroups))}
	
	_nameModes = {"long":0, "medium":1, "short":2} #choices for name lenght
	_defaultNameMode = "short" #short name by default
	
	
	def __init__(self, aminoAcid):
		"""
		Creates an AminoAcid object representing one of the possible Amino Acids.
		aminoAcid can be the name of an amino acid, or an AminoAcid object (in which case a copy is created).
		"""
		self._id = None	#id of the amino acid within the name group
		
		if isinstance(aminoAcid, str):
			self._id = self.__getIdByName(aminoAcid) #id found from aminoAcid name
		elif isinstance(aminoAcid, AminoAcid):
			self._id = aminoAcid._id #copy of id
		else:			
			raise TypeError("aminoAcid must be a string or an AminoAcid object")
		
	
	@staticmethod
	def __getIdByName(name):
		try:
			return AminoAcid._nameDict[name] #get index of name mode
		except:
			raise ValueError("Could not find amino acid name {}".format(name))
			
	@staticmethod
	def getAllNames(nameMode=_defaultNameMode):
		try:
			nameIndex = AminoAcid._nameModes[nameMode] #get index of name mode
		except:
			raise TypeError("nameMode must be 'short', 'medium' or 'long'")
			
		for aa in AminoAcid._nameGroups[1:]: #we exclude the gap (first item)
			yield aa[nameIndex]
	
	#Representation
	def __repr__(self):
		return self._nameGroups[self._id][self._defaultNameMode] #default name mode
		
	def __str__(self):
		return self.getName() #default name mode
	
	def getName(self, nameMode=_defaultNameMode):
		try:
			nameIndex = AminoAcid._nameModes[nameMode] #get index of name mode
		except:
			raise TypeError("nameMode must be 'short', 'medium' or 'long'")
		
		return self._nameGroups[self._id][nameIndex]
	
	#Comparison		
	def __eq__(self, other):
		return self._id == other._id
	
	def __ne__(self, other):
		return self._id != other._id
	
	
	#Hashing
	def __hash__(self):
		return hash(self._id)
