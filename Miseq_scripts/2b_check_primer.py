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
	j = 0; i = 0
	out = open(seqfile.split('.fas')[0]+'_primer.fas','w+')
	for Seq in SeqIO.parse(seqfile,'fasta'):
		K = 0; i += 1; Primerlist = []
		for Primer in SeqIO.parse('SAR_db/Primer_Sequences.fas','fasta'):
			while str(Primer.seq) in str(Seq.seq) and Primer.id not in Primerlist:
				print(Primer.id, " == " , Seq.id)
				Primerlist.append(Primer.id)
				K += 1
# 				for Primer2 in SeqIO.parse('SAR_db/Primer_Sequences.fas','fasta'):
# 					if str(Primer2.seq) in str(Seq.seq) and Primer2.id not in Primerlist:
# 						print(Primer2.id, " =/= " , Seq.id)
# 						K =+ 2

		if K == 2:
			out.write('>'+Seq.description + '\n' + str(Seq.seq) + '\n')
		elif K == 1 :
			j += 1
			print(Seq.id, " matches only one primer")
		elif K == 0:
			j += 1
			print(Seq.id, " does not match the primer")
		else:
			j += 1
			print(Seq.id, " contains more than the 2 primers")
	print(str(j), " sequences do not have both primers of ", str(i), " sequences")			
		
	
def main():
	script, seqfile = argv
	countread(seqfile)
main()