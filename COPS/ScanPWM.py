from Region import Occurence
import re
class ScanPWM:
	@staticmethod
	def scanPwms(nuc,motifs,chro="X",start=0,mm_order=1,mo_min_len=5):
		occurs=[]
		for i in range(mm_order,len(nuc)-mm_order):
			for motif in motifs:
				k=motif.len
				if re.search("N",nuc[i:i+k]) or len(nuc[i:i+k])<k:
					break
				score_p=motif.getScore(nuc[i:i+k],True)
				if score_p > motif.thres_sc:
					occurence=Occurence(chro,i+start,i+k+start,nuc[i-mm_order:i+k],
					score_p,motif.name)
					occurs.append(occurence)
				
				rev_nuc=reverseCompliment(nuc[i:i+k])
				score_n=motif.getScore(rev_nuc,True)
				if score_n > motif.thres_sc:
					occurence=Occurence(chro,i+start,i+k+start,reverseCompliment(nuc[i:i+k+mm_order]),score_n,motif.name)
					occurs.append(occurence)

		return occurs
def reverseCompliment(dna):
	rev_dna=dna[::-1]
	rev_dna=rev_dna.replace('A','1').replace('C','2').replace('T','A').replace('1','T').replace('G','C').replace('2','G')		
	return rev_dna
