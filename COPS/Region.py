class Region:
	def __init__(self,chro=None,begin=None,end=None):
		if chro!=None:		
			self.chr=str(chro)
		else:
			self.chr=""

		if begin!=None:
			self.begin=int(begin)
		else:
			self.begin=0

		if end!=None:
			self.end=int(end)
		else:
			self.end=-1
		
		self.type=None
		self.strand=None
		self.name=None
		self.transcript=None
		self.count=0
		self.reads=[]
		self.source=None

##toString
	def toKey(self):
		if type(self.chr).__name__=='str':
			return self.chr+":"+str(self.begin)+".."+str(self.end)
		else:
			return str(self.chr)+":"+str(self.begin)+".."+str(self.end)
	def toString(self):
		if type(self.chr).__name__=='str':
			return self.chr+"\t"+str(self.begin)+"\t"+str(self.end)
		else:
			return str(self.chr)+"\t"+str(self.begin)+"\t"+str(self.end)
	def toWig(self,read_len):
		##
		counts=[]
		if len(self.reads)>0:
			wigs={}
			for read in self.reads:
				for pos in range(read,read+read_len):
					wigs[pos]=wigs.get(pos,0)+1 
			for pos in range(self.begin,self.end):
				counts.append(wigs.get(pos,0))
		else:
			counts=[ 0 for pos in range(self.begin,self.end)]
	
		return counts

##functional methods
	## check if the reads is in the Region
	def isIn(self,read_pos,read_len=0,perc=0):
		## if 60% of the reads fall into the region, then TRUE
		if int(read_pos+read_len*perc)<self.begin:
			return False
		else:
			if int(read_pos-perc*read_len)<=self.end:
				return True
			else:
				return False
	## check if the paired reads is in the Region
	def isIn_pair(self,read_pos,perc=1):
		## if paired read center fall into the region, then TRUE
		read_len=(read_pos[1]-read_pos[0])
		if int(read_pos[0]+perc*read_len)<self.begin:
			return False
		else:
			if int(read_pos[0]+perc*read_len)<=self.end:
				return True
			else:
				return False

	def isOverlap(self,other):
		##suppose same chromosome
		if(other.begin>=self.begin and other.begin<=self.end):
			return True
		elif(other.end>=self.begin and other.end<=self.end):
			return True
		elif(other.begin<=self.begin and other.end>=self.end):
			return True
		else:
			return False
	def dist(self,other):
		##suppose same chromosome
		##calculate overlap region
		if(other.begin>=self.begin and other.begin<=self.end):
			if other.end > self.end:
				return self.end-other.begin
			else:
				return other.end-other.begin
		elif(other.end>=self.begin and other.end<=self.end):
			if other.begin < self.begin:
				return other.end-self.begin
			else:
				return other.end-other.begin
		elif(other.begin<=self.begin and other.end>=self.end):
			return other.end-other.begin
		else:
			return -1
	def dist_percent(self,other):
		##suppose same chromosome
		## calculate overlap percentages
		if(other.begin>=self.begin and other.begin<=self.end):
			if other.end > self.end:
				return (self.end-other.begin+1)*1.0/self.getLength()
			else:
				return (other.end-other.begin+1)*1.0/self.getLength()
		elif(other.end>=self.begin and other.end<=self.end):
			if other.begin < self.begin:
				return (other.end-self.begin+1)*1.0/self.getLength()
			else:
				return (other.end-other.begin+1)*1.0/self.getLength()
		elif(other.begin<=self.begin and other.end>=self.end):
			return (self.end-self.begin+1)*1.0/other.getLength()
		else:
			return 0
	def dist_2_other(self,other):
		## suppose same chromosome
		## calculate distance between tss
		dist=None
		if self.strand=="+":
			if other.strand=="+":
				dist=self.begin-other.begin
			else:
				dist=self.begin-other.end
		else:
			if other.strand=="+":
				dist=self.end-other.begin
			else:
				dist=self.end-other.end
		return dist

	def getLength(self):
		return int(self.end)-int(self.begin)+1
	def upstream_of_gene(self,gene,extend_up,extend_down):
		result=False
		if gene.strand=="+":
			pbegin=gene.begin-extend_up
			pend=gene.begin+extend_down
			if (pend>gene.end):
				pend=gene.end
			if self.isOverlap(Region(gene.chr,pbegin,pend)):
				result=True
		else:
			pbegin=gene.end-extend_down
			if pbegin<gene.begin:
				pbegin=gene.begin		
			pend=gene.end+extend_up
			if self.isOverlap(Region(gene.chr,pbegin,pend)):
				result=True
		return result
		

class Sequence(Region):
	def __init__(self,chro=None,begin=None,end=None,nuc=None):
		Region.__init__(self,chro,begin,end)
		if nuc!=None:
			self.nuc=nuc
		else:
			self.nuc=""
		self.upstream=None
		self.downstream=None
		self.intron=[]

class Occurence(Region):
	def __init__(self,chr=None,begin=None,end=None,dna=None,score=None,name=None):
		Region.__init__(self,chr,begin,end)
		
		if dna!=None:
			self.dna=dna
		else:
			self.dna=""
		
		if score!=None:
			self.score=score
		else:
			self.score=0
		
		if name!=None:
			self.name=name
		else:
			self.name=""



