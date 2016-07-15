import re
from PWM import PWM
from Region import *

class COPS_IO:
	@staticmethod
	def readMatrix(infile,h_or_v="h",sep="\t"):
		# if h_or_v="h" the input matrix is horizontal, if "v" then it is vertical
		if h_or_v=="h":
			#print self.source
			return getHorizontalMatrix(infile,sep)
		elif h_or_v=="v":
			return getVerticalMatrix(infile,sep)
	@staticmethod
	def readSequences(infile):
		f=open(infile,'r');
		
		seqs=[]
		pattern=re.compile(r'loc=(.*):(\d+)..(\d+)');
		seq=Sequence()
		nuc=""

		for line in f.readlines():
			line=re.sub('\n','',line);
			if re.search("^>",line):
				if seq.chr!=None and len(nuc)>0:
					seq.nuc=nuc.upper()
					if seq.begin==0:
						seq.end=len(nuc)
					nuc=""
					seqs.append(seq)
					seq=Sequence()
				
				se_results=pattern.search(line)
				if se_results!=None:
					results=se_results.groups()
					if len(results)==3:
						seq.chr=results[0]
						seq.begin=int(results[1])
						seq.end=int(results[2])
				else:
					seq.chr="Random"
					seq.begin=0
			else:
				nuc+=line

		seq.nuc=nuc.upper()
		if seq.begin==0:
			seq.end=len(nuc)
		seqs.append(seq)

				
		return seqs
##read horizontal matrix
# ACGT
#>name=tmp
#10	10	10	10	10
#0	0	0	0	0
#0	0	0	0	0
#0	0	0	0	0
def getHorizontalMatrix(source,sep):
	f=open(source,"r")
	pattern=re.compile(r'name=(.*)');
		
	motifs=[]
	motif=PWM()
	matrix=[]

	for line in f.readlines():
		line=re.sub('\n','',line);
		if re.search("^>",line):
			if motif.name!=None and len(matrix)>0:
				motif.setMatrix(matrix)
				matrix=[]
				motifs.append(motif)
				motif=PWM()
			
			results=pattern.search(line).groups()
			if len(results)==1:		
				motif.name=results[0]
		elif re.search("\d+",line):
			tmp=line.split(sep)
			if len(matrix)==0:
				for cc in tmp:
					matrix.append([float(cc)])
			else:
				for i in range(len(matrix)):
					matrix[i].append(float(tmp[i]))
	
	motif.setMatrix(matrix)
	motifs.append(motif)
	return motifs

## read vertical matrix
#ACGT
#>name=tmp
#0	0	0	0
#0	0	0	0
#0	0	0	0
#0	0	0	0
#0	0	0	0
def getVerticalMatrix(source,sep):
	f=open(source,"r")
	pattern=re.compile(r'name=(.*)');
	
	motifs=[]
	motif=PWM()
	matrix=[[],[],[],[]]

	for line in f.readlines():
		line=re.sub('\n','',line);
		if re.search("^>",line):
			if motif.name!=None and len(matrix)>0:
				motif.setMatrix(matrix)
				matrix=[[],[],[],[]]
				motifs.append(motif)
				motif=PWM()
				
			results=pattern.search(line).groups()
			if len(results)==1:		
				motif.name=results[0]
		elif re.search("\d+",line):
			print line
			tmps=line.split(sep)
			for index in range(len(tmps)):
				print tmps[index]
				matrix[index].append(tmps[index])
				motif.setMatrix(matrix)
	motifs.append(motif)
	return motifs
