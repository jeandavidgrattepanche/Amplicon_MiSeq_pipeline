#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "3, March 2, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys, os, re, time, string, os.path
from distutils import spawn
from Bio import SeqIO
from sys import argv

vsearch_path = spawn.find_executable("vsearch")
SAR_db = "Amoeba_db/Am_Cleaned3_5_25_clean.fasta" ##change to the correct database
blastdict = {}

def getBLAST( NGSfile, idmin, qcov, Taxa, readcutoff):
	print("start BLAST SSU_Am_plusout")
	ublast_self = vsearch_path + ' --usearch_global '+NGSfile+' --db '+SAR_db+ ' --strand plus --id '+str(idmin/100)+' --query_cov '+ str(qcov/100)+' --blast6out temp/output.tsv ' ## No -evalue 1e-15 as usearch
	os.system(ublast_self)
	for blast_record in open('temp/output.tsv','r'):
		blastdict[blast_record.split('\t')[0]] = blast_record.split('\n')[0]

	outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch.fasta','w+')
	outseqAm = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_Amoeba_vsearch.fasta','w+')
	for seq in SeqIO.parse(NGSfile,'fasta'):
		try:
			blastdict[seq.id]
			ident =  blastdict[seq.id].split('\t')[2]
			seqmatch = float(blastdict[seq.id].split('\t')[3]) -float(blastdict[seq.id].split('\t')[5])
			cnt = blastdict[seq.id].split('\t')[3] 
			Sim = round(float(seqmatch-int(blastdict[seq.id].split('\t')[4])) / float(cnt) * 100)
			ID = blastdict[seq.id].split('\t')[1]
			seqused = 1 + int(blastdict[seq.id].split('\t')[7]) - int(blastdict[seq.id].split('\t')[6])
			cov = round(float(seqused) / float(len(seq.seq)) * 100)
			outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch.fasta','a')
			if int(seq.description.split('_')[1].replace('r','')) > (int(readcutoff)-1):
				outseq.write('>'+seq.description+ '_'+ ID.split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
				outseq.close()

			if ID.split('_')[0] == str(Taxa):# or ID.split('_')[1] == Taxa:
				if int(seq.description.split('_')[1].replace('r','')) > (int(readcutoff)-1):
					print(seq.id, 'blasted with', ID.split(';size=')[0] , " at ", ident , "% and coverage:", cov )
					outseqAm = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_Amoeba_vsearch.fasta','a')
					outseqAm.write('>'+seq.description+ '_'+ ID.split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
					outseqAm.close()
# 				else:
# 					print(seq.id, 'blasted with', ID.split(';size=')[0] , " at ", ident , "% and coverage:", cov , " BUT low read ", int(seq.description.split('_')[1].replace('r','')))
# 					outseqAml = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_Amoeba_vsearch_lread.fasta','a')
# 					outseqAml.write('>'+seq.description+ '_'+ ID.split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
# 					outseqAml.close()
# 			else:
# 				print(seq.id, 'blasted with', ID.split(';size=')[0] , " at ", ident , "% and coverage:", cov , " BUT not ", Taxa)
# 				outseqnAm = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_Amoeba_vsearch_notAm.fasta','a')
# 				outseqnAm.write('>'+seq.description+ '_'+ ID.split('_rid_')[0] + '_' +str(cov)+'_'+ str(Sim) + '%\n'+str(seq.seq) + '\n')
# 				outseqnAm.close()
				

		except:
			print("NO BLAST for ",seq.id)
			outseq = open('outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch.fasta','a')
			if int(seq.description.split('_')[1].replace('r','')) > (int(readcutoff)-1):
				outseq.write('>'+seq.description+ '_No_BLASTrecord\n'+str(seq.seq) + '\n')
				outseq.close()
		
def main():
	script,  NGSfile, idminy, qcovz, Taxa, readcutoff = argv
	getBLAST(NGSfile, float(idminy),float(qcovz), Taxa, readcutoff)
main()
