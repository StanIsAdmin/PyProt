from copy import deepcopy as deepcopy

nameGroups = (
		("none/gap", "gap", "|"),
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
	Represents one of the twenty amino acids
	"""
	#Tuple of tuples listing amino acid names
	_nameGroups = deepcopy(nameGroups)
	
	#Dictionary mapping name to id
	_nameDict = {nameGroups[id][i]:id for i in range(3) for id in range(len(nameGroups))}
	
	_nameModes = {"long":0, "medium":1, "short":2} #choices for name lenght
	_defaultNameMode = "short" #short name by default
	
	
	
	
	def __init__(self, name):
		"""
		Creates an AminoAcid object representing the amino acid with the name 'name'.
		Name can be the full name, or the three or one-letter name of the amino acid.
		"""
				
		if not isinstance(name, str):
			raise TypeError("name must be a string")
		
		self._id = self.__getIdByName(name) #id of the amino acid (private)
	
	@staticmethod
	def __getIdByName(name):
		try:
			return AminoAcid._nameDict[name] #get index of name mode
		except:
			raise ValueError("Could not find amino acid name {name}".format(name=name))
			
	@staticmethod
	def getAllNames(nameMode=_defaultNameMode):
		try:
			nameMode = AminoAcid._nameModes[nameMode] #get index of name mode
		except:
			raise TypeError("nameMode must be 'short', 'medium' or 'long'")
			
		for aa in AminoAcid._nameGroups[1:]: #we exclude the gap (first item)
			yield aa[nameMode]
	
	#Representation
	def __repr__(self):
		return self._nameGroups[self._id][self._defaultNameMode] #default name mode
		
	def __str__(self):
		return self.getName() #default name mode
	
	def getName(self, nameMode=_defaultNameMode):
		try:
			nameMode = AminoAcid._nameModes[nameMode] #get index of name mode
		except:
			raise TypeError("nameMode must be 'short', 'medium' or 'long'")
		
		return self._nameGroups[self._id][nameMode]
	
	#Comparison	
	def __gt__(self, other):
		if not isinstance(other, AminoAcid):
			raise TypeError("Can not compare AminoAcid with other types")
		return self._id > other._id
	
	def __lt__(self, other):
		if not isinstance(other, AminoAcid):
			raise TypeError("Can not compare AminoAcid with other types")
		return self._id < other._id
	
	def __eq__(self, other):
		return not (self < other or self > other)
	
	def __ne__(self, other):
		return self < other or self > other
	
	def __ge__(self, other):
		return not self < other
	
	def __le__(self, other):
		return not self > other
	
	
	#Hashing
	def __hash__(self):
		return hash(self._id)
		
"""
"""