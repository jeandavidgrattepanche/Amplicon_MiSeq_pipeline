#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "1, June 26, 2018"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
import string
import os.path
from math import log10
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'jgrattepanche@smith.edu'
from sys import argv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

in_db = "Amoeba_db/onlyAm_ref_072618.fasta"
out_db = "Amoeba_db/out_db_BLAST.fasta"
Amblastdict = {}; outblastdict = {}

def getBLAST( NGSfile, idmin, Emin, readcut, cutoff):
	print("start BLAST SAR_db")
	BLASTtable = open('outputs/taxonomic_assignment/BLAST_table.tsv','w')
	BLASTtable.write('OTU\treadnumber\tin_ID\tin_similarity\tin_Evalue\tin_coverage\tout_ID\tout_similarity\tout_Evalue\tout_coverage\tEvalueratio(am/out)\tAssessment\n')
	outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted_diff.fasta','w')
	for seq in SeqIO.parse(NGSfile,'fasta') :
		if int(seq.id.split('_')[1].replace('r','')) > readcut:
			seqtemp = open("temp/seq.fas",'w')
			seqtemp.write(seq.format("fasta"))
			seqtemp.close()
			inblastn_cline = NcbiblastnCommandline(task="blastn", query="temp/seq.fas", db=in_db, outfmt="5", evalue = Emin, perc_identity=idmin, num_threads =2, max_hsps=1, out= "temp/in_BLAST.xml")
		#		blast_records = NCBIXML.parse(result_handle)
			stdout, stderr = inblastn_cline()
			inblast_records = NCBIXML.parse(open("temp/in_BLAST.xml"))
			outblastn_cline = NcbiblastnCommandline(task="blastn", query="temp/seq.fas", db=out_db, outfmt="5", evalue = Emin, perc_identity=idmin, num_threads =2, max_hsps=1, out= "temp/out_BLAST.xml")
		#		blast_records = NCBIXML.parse(result_handle)
			stdout, stderr = inblastn_cline()
			inblast_records = NCBIXML.parse(open("temp/in_BLAST.xml"))
			for inblast_record in inblast_records:
				if inblast_record.descriptions:
					for i in range(1):
		#					print (i)
						inident =  inblast_record.alignments[i].hsps[0].identities
						inmatch = inblast_record.alignments[i].hsps[0].match
						inMatchlist = []
						inseqmatch = len(inblast_record.alignments[i].hsps[0].query)-int(inblast_record.alignments[i].hsps[0].query.count('-'))
						cnt = 0
						for inMatch in inmatch:
							inMatchlist.append('x')
						incnt = inMatchlist.count('x') 
						inevalue = inblast_record.alignments[i].hsps[0].expect
						inSim = round(float(inident) / float(incnt) * 100)
						inID = inblast_record.alignments[i].title
						incov = round(float(inseqmatch) / float(len(seq.seq)) * 100)
						print(seq.id, 'blasted with', inID.split(' ')[1].split('_rid_')[0] , " at ", inSim, "%", "E-value:", inevalue, "coverage:", incov )
		#					outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted.fasta','a')
		#					outseq.write('>'+seq.description.replace(';size=','-').replace(';','r')+ '_'+ ID.split(' ')[1].split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
		#					outseq.close()
				else:
					print("NO BLAST for ", inblast_record.query)
					inevalue = 1
					inSim = "NA"
					inID = "NA NA"
					incov = "NA"
		#				outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted.fasta','a')
		#				outseq.write('>'+seq.description.replace(';size=','-').replace(';','r')+ '_NoSARblast\n'+str(seq.seq) + '\n')
		#				outseq.close()
			stdout, stderr = outblastn_cline()
			outblast_records = NCBIXML.parse(open("temp/out_BLAST.xml"))
			for outblast_record in outblast_records:
				if outblast_record.descriptions:
					for i in range(1):
		#					print (i)
						outident =  outblast_record.alignments[i].hsps[0].identities
						outmatch = outblast_record.alignments[i].hsps[0].match
						outMatchlist = []
						outseqmatch = len(outblast_record.alignments[i].hsps[0].query)-int(outblast_record.alignments[i].hsps[0].query.count('-'))
						outcnt = 0
						for outMatch in outmatch:
							outMatchlist.append('x')
						outcnt = outMatchlist.count('x') 
						outevalue = outblast_record.alignments[i].hsps[0].expect
						outSim = round(float(outident) / float(outcnt) * 100)
						outID = outblast_record.alignments[i].title
						outcov = round(float(outseqmatch) / float(len(seq.seq)) * 100)
						print(seq.id, 'blasted with', outID.split(' ')[1].split('_rid_')[0] , " at ", outSim, "%", "E-value:", outevalue, "coverage:", outcov )
				else:
					print("NO BLAST for ", outblast_record.query)
					outevalue = 1
					outSim = "NA"
					outID = "NA NA"
					outcov = "NA"
#			print(float(float(inevalue)/float(outevalue)))
			if log10(float(inevalue)/float(outevalue)) < float(cutoff) :
				result = "in"
				outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted_diff.fasta','a')
				outseq.write('>'+seq.description+ '_'+ str(int(round(log10(float(inevalue)/float(outevalue)),0))) + '\n'+str(seq.seq) + '\n')
				outseq.close()
			if 	log10(float(inevalue)/float(outevalue)) > float(cutoff) :
				result = "out"
	
			BLASTtable = open('outputs/taxonomic_assignment/BLAST_table.tsv','a')
			BLASTtable.write(seq.id.split('_')[0]+'\t'+seq.id.split('_')[1].replace('r','')+'\t'+str(inID.split(' ')[1])+'\t'+str(inSim)+'\t'+str(inevalue)+'\t'+str(incov)+'\t'+str(outID.split(' ')[1])+'\t'+str(outSim)+'\t'+str(outevalue)+'\t'+str(outcov)+'\t'+str(inevalue/outevalue)+'\t'+str(result)+'\n')
			BLASTtable.close()
# 			outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted_diff.fasta','a')
# 			outseq.write('>'+seq.description+ '_'+ str(result) + '\n'+str(seq.seq) + '\n')
# 			outseq.close()


		
def main():
	script,  NGSfile, idminy, Eminz, readcutz, cutoffz = argv
	idmin = float(idminy)
	Emin = float(Eminz)
	cutoff = float(cutoffz)
	readcut = int(readcutz)
	getBLAST(NGSfile, idmin, Emin, readcut, cutoff)
main()