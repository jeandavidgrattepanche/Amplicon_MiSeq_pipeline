#python3 MiSeq_pipeline_SAR_SWARM_part3.py List_samples.txt RaXML_tree #(dataname) 
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



def removeoutgroup(outputpath, outgrouptree, tree, listsample):#, dataname):
	print("Remove OTUs based on outgroup tree\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): 
		print('ERROR! No outgroup_removal folder. Run part 2!')	
		main()
	os.system('python3 Miseq_scripts/8_remove_outgroup_from_tree.py outputs/outgroup_removal/' +outgrouptree + ' outputs/outgroup_removal/' +tree + ' outputs/OTUs/SWARM_postout.txt outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta')
	os.system('python3 Miseq_scripts/9_randomly_subsample_ingroup.py outputs/outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt '+ listsample)# + ' '+dataname)
	os.system('python3 Miseq_scripts/10_taxonomy_Treev2.py outputs/outgroup_removal/' +tree )
	os.system('python3 Miseq_scripts/11_makeOTUtable_ingroup_v2.py outputs/outgroup_removal/SWARM_postout_nosingleton_nochimeras_in_only.txt outputs/outgroup_removal/subsamples.txt '+ listsample)
	os.system('python3 Miseq_scripts/12_createFinalfiles_v2.2.py outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch.fasta outputs/taxonomic_assignment/taxonomy_by_Tree.txt outputs/OTUs_ingroup/SWARM_postout_nosingleton_nochimeras_in_only_subsampled.txt '+listsample) #+ ' '+ dataname)		
	
def main():
	samplefile = sys.argv[1]
	treefile = sys.argv[2]
	dname = sys.argv[3]

	pathA = os.getcwd()
	path = pathA + "/Rawdata/"

	try:
		listsample = samplefile
	except ValueError:
		samplefile = ""	
	if samplefile == "":
		print ('Your input samplefile is empty.  Try again. ')

	try:
		tree = treefile
		outgrouptree = treefile.split('.tre')[0]+'_outgroup.tre'
	except ValueError:
		treefile = ""	
	if treefile == "":
		print ('Your input treefile is empty.  Try again. ')

#	d = input('where is the name of your run? (Can be found in the name of the sequences files e.g. M00763) \n ')
# 	try:
# 		dataname = dname
# 	except ValueError:
# 		dname = ""	
# 	if dname == "":
# 		print ('Your input is empty.  Try again. ')

	outputpath = pathA + '/outputs/'
	if not os.path.exists(outputpath): 
		print('ERROR! No output folder. run part 2!')
		main()
	removeoutgroup(outputpath, outgrouptree, tree, listsample)#, dataname)
main()
