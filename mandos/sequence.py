import sys,os

class Sequence:
	def __init__(self,name="",seq=""):
		self.name = name
		self.seq = seq
		self.qualstr = ""
		self.qualarr = []
	
	#offset of 33 assumed
	def set_qualstr(self,qual):
		self.qualstr = qual
		if len(self.qualarr) == 0:
			for j in self.qualstr:
				self.qualarr.append(ord(j)-33)
	
	#offset of 33 assumed
	def set_qualarr(self,qual):
		self.qualarr = qual
		if len(self.qualstr) == 0:
			for j in self.qualarr:
				self.qualstr += chr(j+33)
	
	def get_fasta(self):
		retstr = ">"
		retstr += self.name
		retstr += "\n"
		retstr += self.seq
		return retstr
	
	def get_fastq(self):
		retstr = "@"
		retstr += self.name
		retstr += "\n"
		retstr += self.seq
		retstr += "\n+\n"
		retstr += self.qualstr
		return retstr
	

if __name__ == "__main__":
	s = Sequence(name="a",seq="acgt")
	print s.name

	
