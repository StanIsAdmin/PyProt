from copy import deepcopy as deepcopy

class AminoAcid:
	"""
	Represents one of the twenty amino acids
	"""
	_nameGroups = (
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
		("undetermined", "und", "X"),
	)
	"""
	Xle or J means either leucine or isoleucine
	Sec or U refers to selenocysteine, 
	and Pyl or O refers to pyrrolysine.
	"""
	
	_defaultNameMode = 2 #short name by default
	
	def __init__(self, name):
		"""
		Creates an AminoAcid object representing the amino acid with the name 'name'.
		Name can be the full name, or the three or one-letter name of the amino acid.
		"""
				
		if not isinstance(name, str):
			raise TypeError("name must be a string")
		if name=="":
			raise ValueError("name cannot be empty")
		
		self._id = None
		
		if len(name) == 1: #short name
			self._id = self.getIdByName(name.upper(), 2)
		elif len(name) == 3: #medium name
			self._id = self.getIdByName(name.lower(), 1)
		else: #long name
			self._id = self.getIdByName(name.lower(), 0)

			
		if self._id == None:
			raise ValueError("could not find amino acid name '%s'"%name)
			
	@staticmethod
	def getAllNames(nameMode=_defaultNameMode):
		for aa in AminoAcid._nameGroups:
			yield aa[nameMode]
	
	@staticmethod
	def getIdByName(name, nameMode=_defaultNameMode):
		id = 0
		for nameGroup in AminoAcid._nameGroups:
			if name == nameGroup[nameMode]:
				return id
			else:
				id += 1
		
		return None
	
	#Representation
	def __repr__(self):
		return self._nameGroups[self._id][self._defaultNameMode] #default name mode
		
	def __str__(self):
		return self.getName()
	
	def getName(self, nameMode="short"):
		if not isinstance(nameMode, str):
			raise TypeError("nameMode must be a string")
		
		nameMode.lower()
		
		if nameMode == "long":
			return self._nameGroups[self._id][0]
		elif nameMode == "medium":
			return self._nameGroups[self._id][1]
		elif nameMode == "short":
			return self._nameGroups[self._id][2]
		else:
			raise ValueError("accepted modes are 'long', 'medium' and 'short'")
			
	def getId(self):
		return self._id
	
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