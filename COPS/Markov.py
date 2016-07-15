import math
import random
import itertools

class Markov:
	def __init__(self,order=1,strand=1,llist=["A","C","G","T"],seudo=1e-20):
		## default order=1
		self.order=order
		self.list=llist
		self.strand=strand
		self.trans_count={}
		self.trans_total={}
		self.trans_prob={}
		if self.order==0:
			for alph in self.list:
				self.trans_count[alph]=0
				self.trans_prob[alph]=seudo
		else:
			for letter in itertools.product(self.list,repeat=order):
				self.trans_count["".join(letter)]={}
				self.trans_total["".join(letter)]=0
				self.trans_prob["".join(letter)]={}
				for alph in self.list:
					self.trans_count["".join(letter)][alph]=0
					self.trans_prob["".join(letter)][alph]=seudo
	def count_trans(self,dna):
		#print dna
		if self.order==0:
			for ind in range(len(dna)):
				if dna[ind] not in self.list:
					continue
				self.trans_count[dna[ind]]+=1
			return 0
		for ind in range(len(dna)-self.order):
			found=True
			for letter in dna[ind:ind+self.order+1]:
				if letter not in self.list:
					found=False
			if found:
				self.trans_count[dna[ind:ind+self.order]][dna[ind+self.order]]+=1
				self.trans_total[dna[ind:ind+self.order]]+=1
	def norm_trans(self,seudo=1e-20):
		if self.order==0:
			for alph in self.list:
				self.trans_prob[alph]=1.0*self.trans_count[alph]/sum(self.trans_count.values())
			return 0
		for key in self.trans_count:
			for alph in self.list:
				if self.trans_count[key][alph]!=0 and self.trans_total[key]!=0:
					self.trans_prob[key][alph]=1.0*self.trans_count[key][alph]/self.trans_total[key]
	def prob_log(self,dna):
		val_log=[]
		if self.order==0:
			for ind in range(len(dna)):
				val_log.append(math.log(self.trans_prob[dna[ind]],2))
			return sum(val_log)
		for ind in range(len(dna)-self.order):
			val_log.append(math.log(self.trans_prob[dna[ind:ind+self.order]][dna[ind+self.order]],2))
		return sum(val_log)
	def toString(self):
		output=[]
		if self.order==0:
			for alph in self.list:
				output.append("\t".join(["0",str(self.strand),alph,str(self.trans_prob[alph])]))
		else:
			for key in self.trans_prob.keys():
				for alph in self.list:
					output.append("\t".join([str(self.order),str(self.strand),key,alph,str(self.trans_prob[key][alph])]))
		return output
class MarkovOB_IO:
	@staticmethod
	def read(infile,order=3,strand=1,llist=["A","C","G","T"]):
		mbg=MarkovOB(order,strand,llist)
		fin=open(infile,"r")
		line=fin.readline()
		while line:
			tmp=line.replace("\n","").split("\t")
			if tmp[0]==str(0):
				mbg.markovs[0].trans_prob[tmp[2]]=float(tmp[3])
			else:
				mbg.markovs[int(tmp[0])].trans_prob[tmp[2]][tmp[3]]=float(tmp[4])
			line=fin.readline()
		return mbg

class MarkovOB:
	def __init__(self,order=3,strand=1,llist=["A","C","G","T"],seudo=1e-20):
		self.markovs=[]
		self.order=order
		self.strand=1
		for i in range(self.order+1):
			mm=Markov(i,strand,llist,seudo)
			self.markovs.append(mm)
	def train(self,nuc):
		for i in range(self.order+1):
			self.markovs[i].count_trans(nuc)
	def read_fasta(self,file_in):
		fin=open(file_in,"r")
		line=fin.readline()
		count=1
		dna=None
		while line:
			if count%1000==0:
				print "reading",count,"lines"
			if line.find(">")!=-1:
				if dna:
					nuc="".join(dna)
					self.train(nuc)
					if self.strand==2:
						self.train(rev_comp(nuc))
				dna=[]
				line=fin.readline()
				count+=1
				continue
			#print line
			dna.append(line.replace("\n","").upper())
			#self.count_trans(line.replace("\n","").upper())
			line=fin.readline()
			count+=1
		self.train("".join(dna))
		if self.strand==2:
			self.train(rev_comp("".join(dna)))
	def norm(self):
		for i in range(self.order+1):
			self.markovs[i].norm_trans()
	def prob_log(self,dna):
		probs=[]
		probs.append(math.log(self.markovs[0].trans_prob[dna[0]],2))
		for i in range(1,self.order):
			probs.append(math.log(self.markovs[i].trans_prob[dna[0:i]][dna[i]],2))
		for i in range(len(dna)-self.order):
			probs.append(math.log(self.markovs[self.order].trans_prob[dna[i:i+self.order]][dna[i+self.order]],2))
		return sum(probs)
	def dist(self,dna,mbg):
		return (self.prob_log(dna)-mbg.prob_log(dna))/len(dna)
	def write_out(self,output):
		fout=open(output,"w")
		for i in range(self.order+1):
			for line in self.markovs[i].toString():
				fout.write(line+"\n")
class MarkovBG:
	def __init__(self,order=1,strand=1,llist=["A","C","G","T"],seudo=1e-20):
		## default order=1
		self.order=order
		self.list=llist
		self.strand=strand
		self.count_kmer={}
		self.kmer_total=0
		if self.order==0:
			for alph in self.list:
				self.count_kmer[alph]=0
		else:
			for letter in itertools.product(self.list,repeat=self.order):
				for alph in self.list:
					self.count_kmer["".join(letter)+alph]=0

	def count(self,dna):
		#print dna
		for ind in range(len(dna)-self.order):
			found=True
			for letter in dna[ind:ind+self.order+1]:
				if letter not in self.list:
					found=False
			if found:
				self.count_kmer[dna[ind:ind+self.order+1]]+=1
				self.kmer_total+=1

	def read_fasta(self,file_in):
		fin=open(file_in,"r")
		line=fin.readline()
		count=1
		dna=None
		while line:
			if count%1000==0:
				print "reading",count,"lines"
			if line.find(">")!=-1:
				if dna:
					nuc="".join(dna)
					self.count(nuc)
					if self.strand==2:
						self.count(rev_comp(nuc))
				dna=[]
				line=fin.readline()
				count+=1
				continue
			#print line
			dna.append(line.replace("\n","").upper())
			#self.count_trans(line.replace("\n","").upper())
			line=fin.readline()
			count+=1
		self.count("".join(dna))
		if self.strand==2:
			self.count(rev_comp("".join(dna)))
	def write_log_bg(self,output):
		fout=open(output,"w")
		log_total=math.log(self.kmer_total,2)
		count_total=len(self.count_kmer.keys())
		sorted_counts=sorted(self.count_kmer.values(),reverse=True)
		p_count=0.0
		for count_num in sorted_counts:
			if (p_count+1)/count_total>0.05:
				break
			fout.write(str(self.order+1)+"\t"+str((math.log(count_num,2)-log_total)/(self.order+1))+"\t"+str(p_count/count_total)+"\n")
			p_count+=1

	def write_all(self,output):
		fout=open(output,"w")
		for key in self.count_kmer.keys():
			fout.write(key+"\t"+str(self.count_kmer[key]*1.0/self.kmer_total)+"\n")

def rev_comp(dna):
	rev_dna=dna[::-1]
	rev_dna=rev_dna.replace('A','1').replace('C','2').replace('T','A').replace('1','T').replace('G','C').replace('2','G')
		
	return rev_dna
			
