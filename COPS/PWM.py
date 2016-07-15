import math
import random
import re
class PWM:
	def __init__(self,name=None,matrix=None):
		if name!=None:
			self.name=name
		else:
			self.name=None

		if matrix!=None:
			self.matrix=matrix
		else:
			self.matrix=[]
		self.pssm=[]
		self.pfm=[]

		self.len=0
		self.total=0
		self.max=0
		self.min=0
		self.r_maxs=[] ## speed up search, score+rest<score, then go next
		self.thres_sc=0
##setter
	def setMatrix(self,matrix):
		self.matrix=matrix
		self.len=len(self.matrix)
		self.total=sum(self.matrix[0])
		self.pfm=[[] for i in range(self.len)]
		for i in range(self.len):
			self.pfm[i]=[self.matrix[i][j]/self.total for j in range(4)]
		self.pwm2pssm()


##2pssm
	def pwm2pssm(self,bg=[0.298,0.202,0.202,0.298],seudo=1,m_type=1):
		## transform counting matrix to pssm
		##correcting probability by adding seudo count

		pssm=[[] for i in range(self.len)]
		maxs=[]
		for i in range(self.len):
			if m_type==1:## use counting matrix
				pssm[i]=[math.log(((self.matrix[i][j]+seudo)/(self.total+4*seudo))/bg[j]) for j in range(4)]
			else:
				pssm[i]=[math.log(((self.pfm[i][j]+seudo)/(self.total+4*seudo))/bg[j],2) for j in range(4)]
			maxs.append(max(pssm[i]))
			self.max+=max(pssm[i])
			self.min+=min(pssm[i])
		self.pssm=pssm
		for i in range(1,len(maxs)):
			self.r_maxs.append(sum(maxs[i:]))
		self.r_maxs.append(0)
		

	def getScore(self,nuc,thres=False):
		index={'A':0,'C':1,'G':2,'T':3}
		score=0
		if not thres:
			for j in range(len(nuc)):
				score=score+self.pssm[j][index[nuc[j]]]
			return score
		#print self.name,self.thres_sc,nuc
		for j in range(len(nuc)):
			score+=self.pssm[j][index[nuc[j]]]
			if score+self.r_maxs[j]<self.thres_sc:
				#print j,score,self.r_maxs[j]
				return 0
		return score
	def calThres(self,p_thres=0.001,bg=[0.298,0.202,0.202,0.298]):
		## init distribution
		dist={}
		dist[0]={}
		for j in range(4):
			dist[0][self.pssm[0][j]]=dist[0].get(self.pssm[0][j],0)+bg[j]
		
		## start calc distribution
		for i in range(1,self.len):
			dist[i]={}
			if i-2>=0:
				## release some memory
				del dist[i-2]
			for pre_sc in dist[i-1].keys():
				for j in range(4):
					sc=pre_sc+self.pssm[i][j]
					dist[i][sc]=dist[i].get(sc,0)+dist[i-1][pre_sc]*bg[j]
		pval=0
		## release memory
		del dist[self.len-2]

		s_sorted=sorted(dist[self.len-1].keys(),reverse=True)
		#s_sorted=sorted(dist[self.len-1].keys())
		for sc in s_sorted:
			#print sc,dist[self.len-1][sc]
			if pval+dist[self.len-1][sc]>p_thres:
			#if pval>(1-p_thres):
				self.thres_sc=sc
				return sc
			pval+=dist[self.len-1][sc]
##toString
	def matrix2string(self,m_type=1):
		results=[]
		for i in range(self.len):
			if m_type==1: ## original count matrix
				results.append("\t".join([str(cc) for cc in self.matrix[i]]))
			elif m_type==2: ## frequenty matrix
				results.append("\t".join([str(ff) for ff in self.pfm[i]]))
			elif m_type==3:
				results.append("\t".join([str(pp) for pp in self.pssm[i]]))
		return "\n".join(results)



