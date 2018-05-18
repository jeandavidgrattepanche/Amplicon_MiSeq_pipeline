#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "6, April 21, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"


import re,os, sys
from Bio import SeqIO
from Bio import Phylo
from sys import argv


def main():
	script, seqfile, cutoff = argv
	IDlist = []; IDdict={}
	out= open('outputs/chimeras/Seq_reads_'+str(cutoff)+'.fas','w+')
	for Seq in SeqIO.parse(seqfile,'fasta'):
		if 'SWARM' in Seq.id:
			read = Seq.id.split('_')[0].split('-')[1].replace('r','')
			if int(read) >= int(cutoff):
				out.write('>'+Seq.id +'\n'+str(Seq.seq)+'\n')
		else:
			out.write('>'+Seq.id +'\n'+str(Seq.seq)+'\n')
	out.close()		
		
main()


##following command:
## python3 Miseq_scripts/7_remove_column.py outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_5000_BLASTed_TA.fasta fasta SAR_db/SSU_SAR_EUK_v14.3_mask_75.txt
## curl -u $CRA_USER:$PASSWORD -H cipres-appkey:$KEY $URL/job/$CRA_USER -F tool=RAXMLHPC8_REST_XSEDE -F vparam.select_analysis_=fv -F input.infile_=@./outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_5000_BLASTed_TA_masked.fas -F vparam.runtime_=168 -F vparam.dna_gtrcat_=GTRGAMMA -F vparam.invariable_=I -F input.treetop_=@./SAR_db/SSU_SAR_EUK_v14.3_RAxML_constraint_rooted.tre -F metadata.statusEmail=true
