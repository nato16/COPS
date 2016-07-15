from __future__ import generators
class Node:
	def __init__(self,key="",count=1,childs=[],parent=None):
		self.key=key
		self.count=count
		self.childs=[]
		self.parent=parent
		self.str=""

## setter
	def setKey(self,key):
		self.key=key
	def addChild(self,node):
		self.childs.append(node)
	def setChilds(self,childs):
		self.childs=childs
	def setParent(self,node):
		self.parent=node
	
## getter
	def getKey(self):
		return self.key
	def getCount(self):
		return self.count
	def getChilds(self):
		return self.childs
	def getParent(self):
		return self.parent

## recursive 
	def bottomUp(self):
		if self.parent==None:
			return self.toString()
		else:
			return self.toString()+"-->"+self.parent.bottomUp()
	def bottomUpKey(self):
		if self.parent==None:
			return {self.key:self.count}
		else:
			return dict({self.key:self.count}.items()+self.parent.bottomUpKey().items())
## toString()
	def toString(self):
		stri=self.key+":"+str(self.count)+",p:"
		if self.parent!=None:
			return stri+self.parent.key
		else:
			return stri+"None"

def uniq_combinations(items, n):
	if n==0: yield []
	else:
		for i in xrange(len(items)):
			for cc in uniq_combinations(items[i+1:],n-1):
				yield [items[i]]+cc

class Tree:
	def __init__(self,pattern_len=2,motif_in_order=[]):
		self.start=Node("start")
		self.level=0
		self.current=None
		self.orders=motif_in_order
		## save connection of same key
		self.pointer={}
		self.pattern_len=int(pattern_len)
		self.key_2_count={}

## insert new node, and update the count{'a':2,'b':3}
	def updateTree(self,dictl,seq=False):
		level=0
		llist=[]
		if seq:
			## count by seq, norm the count to 1
			for key in dictl.keys():
				dictl[key]=1
		for ol in self.orders:
			if ol in dictl.keys():
				llist.append(ol)
		if len(self.start.childs)==0:
			## empty tree, simply add one by one
			self.current=self.start
			for item in llist:
				level+=1
				tmp=Node(item,dictl[item])
				self.current.addChild(tmp)
				tmp.setParent(self.current)
				## add to the pointer
				self.pointer[item]=[]
				self.pointer[item].append(tmp)
				self.current=tmp
		else:
			## non empty, compare then add
			self.current=self.start
			level=0
			for item in llist:
				level+=1
				found=False
				for child in self.current.getChilds():
					if child.key==item:
						found=True
						child.count+=dictl[item]
						self.current=child
				if found==False:
					tmp=Node(item,dictl[item])
					self.current.addChild(tmp)
					tmp.setParent(self.current)
					## add to the pointers
					if item not in self.pointer.keys():
						self.pointer[item]=[]
					self.pointer[item].append(tmp)
					self.current=tmp
				
		if level>self.level:
			self.level=level
	def walkTree(self):
		for child in self.start.childs:
			print child.key,child.count
			for child2 in child.childs:
				print ">",child2.key,child2.count
				for child3 in child2.childs:
					print ">>",child3.key,child3.count
					for child4 in child3.childs:
						print ">>>",child4.key,child4.count
						for child5 in child4.childs:
							print ">>>>",child5.key,child5.count

## getter
	def getLevel(self):
		return self.level	
	def getStart(self):
		return self.start
	def getPointer(self):
		return self.pointer
## generate patterns
	def generate_pattern(self,keyVals,primary_key):
		#{'a': 7, 'c': 3, 'b': 5,'start':1}
		del keyVals['start']
		#print primary_key
		if len(keyVals.keys())<self.pattern_len:
			return None
		#print primary_key,self.pattern_len,keyVals
		pat_count={}
		
		#print uniq_combinations(sorted(keyVals.keys()),self.pattern_len)
		for pat in uniq_combinations(sorted(keyVals.keys()),self.pattern_len):
			#print pat
			if primary_key in pat:
				#print ";".join(pat)
				pat_count[';'.join(pat)]=min([keyVals[key] for key in pat])
		return pat_count
		#print keyVals,pat_count
	def parse_tree(self,number=1):
		#print "parsing the data"
		results={}
		for item in self.orders:
			if item in self.pointer.keys():
				#print item
				tmp1={}
				for node in self.pointer[item]:
					keys=node.bottomUpKey()
					#if "twi" in keys and "bap_fur" in keys:
					#	print keys["twi"],keys["bap_fur"]
					#print "\t",keys
					tmp2=self.generate_pattern(keys,node.key)
					if tmp2==None:
						continue	
					for key in tmp2.keys():
						if key not in tmp1.keys():
							tmp1[key]=int(tmp2[key])
						else:
							tmp1[key]+=int(tmp2[key])
				for key in tmp1.keys():
					self.key_2_count[key]=tmp1[key]
					if tmp1[key]<number:
						del tmp1[key]
					else:
						results[key]=tmp1[key]
		return results
