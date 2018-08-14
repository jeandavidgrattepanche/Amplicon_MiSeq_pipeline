#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "1, August 13, 2018"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
from Bio import SeqIO
from sys import argv

def countread(seqfile):	
	j = 0; i = 0; bases= ['A','G','T','C']
	out = open(seqfile.split('.fas')[0]+'_primer.fas','w+'); out2 = open(seqfile.split('.fas')[0]+'_NOprimer.fas','w+')
	for Seq in SeqIO.parse(seqfile,'fasta'):
		K = 0; i += 1; Primerlist = []; newseq = str(Seq.seq)
		for Primer in SeqIO.parse('SAR_db/Primer_Sequences.fas','fasta'):
			for l in range(0,len(str(Primer.seq))):
				for b in bases:
					newprim = str(Primer.seq)[0:l] + b + str(Primer.seq)[l+1:]
# 					print(str(Primer.seq), '  ', newprim)
					for m in range(0,len(str(Primer.seq))):
						for d in bases:
							newprim = str(Primer.seq)[0:m] + d + str(Primer.seq)[m+1:]
							if newprim in str(Seq.seq) and Primer.id.split('_')[0] not in Primerlist:
#								print(Primer.id, " == " , Seq.id)
								Primerlist.append(Primer.id.split('_')[0])
								K += 1
								if "f" in Primer.id:
									newseq = newseq.split(str(newprim))[1]
								if "r" in Primer.id:
									newseq = newseq.split(str(newprim))[0]
		if K == 2:
			out = open(seqfile.split('.fas')[0]+'_primer.fas','a')
			out.write('>'+Seq.description + '\n' + newseq + '\n')
			out.close()
		else:
			j += 1
			out2 = open(seqfile.split('.fas')[0]+'_NOprimer.fas','a')
			out2.write('>'+Seq.description + '\n' + str(Seq.seq) + '\n')
			out2.close()
		print(str(j), " sequences do not have both primers of ", str(i), " sequences" , end= '\r')
		
	
def main():
	script, seqfile = argv
	countread(seqfile)
main()