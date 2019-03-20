#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "8, June 1st, 2018"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []



def removeoutgroup(outputpath, outgrouptree, tree, listsample, dataname):
	print("Remove OTUs based on outgroup tree\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): 
		print('ERROR! No outgroup_removal folder. Run part 2!')	
		main()
	os.system('python3 Miseq_scripts/8_remove_outgroup_from_tree.py outputs/outgroup_removal/' +outgrouptree + ' outputs/outgroup_removal/' +tree + ' outputs/OTUs/SWARM_postout.txt outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta')
	os.system('python3 Miseq_scripts/9_randomly_subsample_ingroup.py outputs/outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt '+ listsample)# + ' '+dataname)
	os.system('python3 Miseq_scripts/10_taxonomy_Treev2.py outputs/outgroup_removal/' +tree )
	os.system('python3 Miseq_scripts/11_makeOTUtable_ingroup_v2.py outputs/outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt outputs/outgroup_removal/resubsamples.txt '+ listsample)
	os.system('python3 Miseq_scripts/12_createFinalfiles_v2.2.py outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch.fasta outputs/taxonomic_assignment/taxonomy_by_Tree.txt outputs/OTUs_ingroup/SWARM_postout_nosingleton_nochimeras_in_only_subsampled.txt '+listsample) #+ ' '+ dataname)		
	
def main():
#  	a =sys.argv[1]
	b = sys.argv[1]
	d = sys.argv[2]
	c = sys.argv[3]
#	a = input('where your raw data folder is (should be a folder:  /Users/katzlab33/Documents/MiSeq2016/MiSeq_pipeline ) \n ')
# 	try:
# 		Path = a
# 	except ValueError:
# 		a = ""	
# 	if a == "":
# 		print ('Your input is empty.  Try again. ')
# 	else:
#		pathA = a.split(' ')[0]
	pathA = os.getcwd()
	path = pathA + "/Rawdata/"

#	b = input('where is your sample list file (should be a file:  samplelist.txt: LKM## (tab) samplename ) \n ')
	try:
		listsample = b
	except ValueError:
		b = ""	
	if b == "":
		print ('Your input is empty.  Try again. ')

#	d = input('where is the name of your run? (Can be found in the name of the sequences files e.g. M00763) \n ')
	try:
		dataname = d
	except ValueError:
		d = ""	
	if d == "":
		print ('Your input is empty.  Try again. ')
		
		
#	print("The trees needs to be in the outgroup_removal folder and in newick format\n")
#	c = input('what is the name of your tree (for example:  RAxML_labelledTree_masked_<project>.tre and the same name with _outgroup at the end for the outgroup tree :RAxML_labelledTree_masked_<project>_outgroup.tre ) \n Do Not add the full path or the script will break \n')
	try:
		tree = c
		outgrouptree = c.split('.tre')[0]+'_outgroup.tre'
	except ValueError:
		c = ""	
	if c == "":
		print ('Your input is empty.  Try again. ')

	outputpath = pathA + '/outputs/'
	if not os.path.exists(outputpath): 
		print('ERROR! No output folder. run part 2!')
		main()
	removeoutgroup(outputpath, outgrouptree, tree, listsample, dataname)
main()
