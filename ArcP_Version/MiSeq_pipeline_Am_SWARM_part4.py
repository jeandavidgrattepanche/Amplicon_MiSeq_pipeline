#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "7.1,April 19, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []



def removeoutgroup(outputpath, outgrouptree, tree, listsample):
	print("Remove OTUs based on outgroup tree\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): 
		print('ERROR! No outgroup_removal folder. Run part 2!')	
		main()
	os.system('python3 Miseq_scripts/8_remove_outgroup_from_tree.py outputs/outgroup_removal/' +outgrouptree + ' outputs/outgroup_removal/' +tree + ' outputs/OTUs/SWARM_postout.txt outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta')
	os.system('python3 Miseq_scripts/9_randomly_subsample_ingroup.py outputs/outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt '+ listsample)
	os.system('python3 Miseq_scripts/10_taxonomy_Treev2.py outputs/outgroup_removal/' +tree )
	os.system('python3 Miseq_scripts/11_makeOTUtable_ingroup_v2.py outputs/outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt outputs/outgroup_removal/subsampled.txt '+ listsample)
	os.system('python3 Miseq_scripts/12_createFinalfiles_diff_v3.py outputs/outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.fasta outputs/taxonomic_assignment/BLAST_table.tsv outputs/taxonomic_assignment/taxonomy_by_Tree.txt outputs/OTUs_ingroup/SWARM_postout_nosingleton_nochimeras_in_only_subsampled.txt '+listsample)		
	
def main():
	b = sys.argv[1]
	listsamp = []
	pathA = os.getcwd()
	path = pathA + "/Rawdata/"
	try:
		listsample = b
	except ValueError:
		b = ""	
	if b == "":
		print ('Your input is empty.  Try again. ')
	for samp in open(listsample,'r'):
		if samp.split('\t')[0] not in listsamp:
			listsamp.append(samp.split('\t')[0])
	c = input('what is the name of your tree (for example:  RAxML_labelledTree_masked_<project>.tree and the same name with _outgroup at the end for the outgroup tree :RAxML_labelledTree_masked_<project>_outgroup.tree ) \n Do Not add the full path or the script will break \n')
	try:
		tree = c
		outgrouptree = c.split('.tree')[0]+'_outgroup.tree'
	except ValueError:
		c = ""	
	if c == "":
		print ('Your input is empty.  Try again. ')

	outputpath = pathA + '/outputs/'
	if not os.path.exists(outputpath): 
		print('ERROR! No output folder. run part 2!')
		main()
	removeoutgroup(outputpath, outgrouptree, tree, listsample)
main()
