from aminoAcid import AminoAcid

class Sequence:
	"""
	Represents an amino acid sequence, and accepts all operations that can be done on such a sequence
	"""
	
	def __init__(self, aminoAcids=None):
		"""
		Creates a Sequence object that represents the amino acid sequence contained in aminoAcids.
		aminoAcids can be one of the following :
		- None, meaning the Sequence is empty (default)
		- an AminoAcid object or a string of valid AminoAcid name
		- a list of AminoAcid objects or strings of valid AminoAcid names
		"""
		
		self._nameMode = "short" #the way in which AA names are displayed

		#Format aminoAcids into a list of AminoAcid objects.
		self._aaList = self.__formatAAList(aminoAcids)
		
	@staticmethod
	def __formatAAList(aminoAcids):
		#Sequence object is deep copied
		if isinstance(aminoAcids, Sequence):
			return deepcopy(aminoAcids)
		
		#None becomes an empty Sequence
		elif aminoAcids is None:
			return []
		
		#A string is converted to a list
		elif isinstance(aminoAcids, str):
			if aminoAcids[0:5] == "short": #Multiple Amino Acids in short name mode
				return [AminoAcid(aa) for aa in aminoAcids[5:]]
			else: #A single Amino Acid with any name mode
				return [AminoAcid(aminoAcids)]
		
		#An AminoAcid is put within a list
		elif isinstance(aminoAcids, AminoAcid):
			return [deepcopy(aminoAcids)]
		
		#Only lists are allowed besides the previous
		elif not isinstance(aminoAcids, list):
			raise TypeError("aminoAcids must be a list, an AminoAcid object or a string")
		
		aaCopy = None
		formattedAAList = []
		for aa in aminoAcids:
			if isinstance(aa, AminoAcid):
				aaCopy = deepcopy(aa) #create a copy of the object
			else:
				aaCopy = AminoAcid(aa) #create the object from its name
			formattedAAList.append(aaCopy)
			
		return formattedAAList
	
	#Size
	def __len__(self):
		return len(self._aaList)
	def __gt__(self, other):
		if not isinstance(other, Sequence):
			raise ValueError("cannot compare Sequence object with non-Sequence object")
		return len(self) > len(other)
	def __lt__(self, other):
		if not isinstance(other, Sequence):
			raise ValueError("cannot compare Sequence object with non-Sequence object")
		return len(self) < len(other)
	def __ge__(self, other):
		if not isinstance(other, Sequence):
			raise ValueError("cannot compare Sequence object with non-Sequence object")
		return len(self) >= len(other)
	def __le__(self, other):
		if not isinstance(other, Sequence):
			raise ValueError("cannot compare Sequence object with non-Sequence object")
		return len(self) <= len(other)
	
	#Comparison
	def __eq__(self, other):
		if not isinstance(other, Sequence):
			raise ValueError("cannot compare Sequence object with non-Sequence object")
		return len(self) == len(other) and all(self[i] == other[i] for i in len(self))
	
	def __ne__(self, other):
		if not isinstance(other, Sequence):
			raise ValueError("cannot compare Sequence object with non-Sequence object")
		return not self == other
	
	#Iteration
	def __iter__(self):
		return iter(self._aaList)
	
	#Representation
	def __repr__(self):
		return "-".join([aa.getName(self._nameMode) for aa in self])
	
	def changeNameMode(self, newMode):
		newMode.lower()
		if newMode in ("long", "medium", "short"):
			self._nameMode = newMode
		else:
			raise ValueError("newMode must be 'long', 'medium' or 'short'")
	
	#Manipulation
	def __getitem__(self, key):
		return deepcopy(self._aaList[key])
		
	def __setitem__(self, key, value):
		self._aaList[key] = AminoAcid(value)
		
	def __delitem__(self, key):
		del self._aaList[key]
	
	def insert(self, sequence, index=None):
		if index is None:
			index = len(self)
		#We need a formatted deep copy of sequence
		for aa in self.__formatAAList(sequence):
			self._aaList.insert(index, aa)
			index += 1
		
	def __contains__(self, item):
		return self.contains(item)
	
	def contains(self, sequence, index=0):
		"""
		Checks if 'sequence' is a sub-sequence of self, starting at index. Returns
		- (start, end) indexes of the sub-sequence within self (end index not included)
		- False if sequence is not a sub-sequence of self
		"""
		if not isinstance(index, int) or index < 0 or index >= len(self):
			raise ValueError("index must be a positive integer lesser than len(self)")
		#We won't change sequence, there is no need to deep copy it
		if not isinstance(sequence, Sequence):
			sequence = self.__formatAAList(sequence)
			
		if len(self)-index < len(sequence): #sub-sequences must be smaller or equal in size
			return False
		
		for i in range(index, len(self)-len(sequence)+1):
			for j in range(len(sequence)):
				if self[i+j] != sequence[j]:
					break
			else:
				return i, i+len(sequence)
		return False
		
	def remove(self, sequence, index=0):
		"""
		Removes the first occurrence of 'sequence' in self, starting at index. Returns
		- (start, end) indexes of the sub-sequence within self (end index not included)
		- False if sequence is not a sub-sequence of self
		"""
		if not isinstance(index, int) or index < 0 or index >= len(self):
			raise ValueError("index must be a positive integer lesser than len(self)")
		
		#We won't change sequence, there is no need to deep copy it
		if not isinstance(sequence, Sequence):
			sequence = self.__formatAAList(sequence)
			
		indexes = self.contains(sequence, index)
		if indexes == False:
			raise ValueError("sequence is a sub-sequence of self after index")
		else:
			del self._aaList[indexes[0]:indexes[1]]
			
	def delete(self, start=0, stop=None):
		"""
		Deletes Amino Acids between indexes start (included) and stop (excluded).
		If start is not specified, deletion will begin at index 0.
		If stop is not specified, deletion will stop after one item.
		"""
		if stop is None:
			stop = start + 1
		if stop <= start:
			raise ValueError("stop index must be strictly greater than start index")
		del self._aaList[indexes[0], indexes[1]]
			
		
"""					
l = ["lysine", "G", AminoAcid("arg")]
s = Sequence(l)
print(l)
print(s)
s.changeNameMode("long")
print(s)
s.insert("ALA", 3)
print(s)
s[0] = "ALA"
print(s)
del s[2]
print(s)
s.remove(["glycine", "ala"])
print(s)
"""