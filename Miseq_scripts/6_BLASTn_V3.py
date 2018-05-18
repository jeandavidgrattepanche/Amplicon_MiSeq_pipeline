#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "3, March 2, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
import string
import os.path
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'jgrattepanche@smith.edu'
from sys import argv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

SAR_db = "SAR_db/SAR_db"

def getBLAST( NGSfile, idmin, Emin):
	print("start BLAST SAR_db")
	outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted.fasta','w')
	for seq in SeqIO.parse(NGSfile,'fasta'):
		seqtemp = open("temp/seq.fas",'w')
		seqtemp.write(seq.format("fasta"))
		seqtemp.close()
		blastn_cline = NcbiblastnCommandline(task="blastn", query="temp/seq.fas", db=SAR_db, outfmt="5", evalue = Emin, perc_identity=idmin, num_threads =4, max_hsps=1, out= "temp/output.xml")
#		blast_records = NCBIXML.parse(result_handle)
		stdout, stderr = blastn_cline()
		blast_records = NCBIXML.parse(open("temp/output.xml"))
		for blast_record in blast_records:
			if blast_record.descriptions:
				for i in range(1):
#					print (i)
					ident =  blast_record.alignments[i].hsps[0].identities
					match = blast_record.alignments[i].hsps[0].match
					Matchlist = []
					seqmatch = len(blast_record.alignments[i].hsps[0].query)-int(blast_record.alignments[i].hsps[0].query.count('-'))
					cnt = 0
					for Match in match:
						Matchlist.append('x')
					cnt = Matchlist.count('x') 
					evalue = blast_record.alignments[i].hsps[0].expect
					Sim = round(float(ident) / float(cnt) * 100)
					ID = blast_record.alignments[i].title
					cov = round(float(seqmatch) / float(len(seq.seq)) * 100)
					print(seq.id, 'blasted with', ID.split(' ')[1].split('_rid_')[0] , " at ", Sim, "%", "E-value:", evalue, "coverage:", cov )
					outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted.fasta','a')
					outseq.write('>'+seq.description.replace(';size=','-').replace(';','r')+ '_'+ ID.split(' ')[1].split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
					outseq.close()
			else:
				print("NO BLAST for ", blast_record.query)
				outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted.fasta','a')
				outseq.write('>'+seq.description.replace(';size=','-').replace(';','r')+ '_NoSARblast\n'+str(seq.seq) + '\n')
				outseq.close()
		
def main():
	script,  NGSfile, idminy, Eminz = argv
	idmin = float(idminy)
	Emin = float(Eminz)
	getBLAST(NGSfile, idmin, Emin)
main()