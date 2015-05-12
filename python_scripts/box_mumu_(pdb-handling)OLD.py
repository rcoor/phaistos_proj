import os
import subprocess
import numpy as np
'Importing the BioPDB package:'
from Bio.PDB import *

def pdb_handling(prot_name,show_pp=0,check_residueNO=False):
	#at the moment this script handles a single chain in multiple chain homomeric proteins


	'file must be available in folder'
	pdb_file = prot_name+".pdb"

	'Uncomment this section if the PDB is to be downloaded'
	# 'PDB download and saving'
	# pdbl = PDBList()
	# pdb_file = pdbl.retrieve_pdb_file(prot_name)

	'.pdb structure is parsed'
	parser = PDBParser()
	structure = parser.get_structure(prot_name, pdb_file)
	
	'prepares a folder for modified .pdbs'
	if not os.path.exists('mod_pdb'):
		os.makedirs('mod_pdb')

	'PDBIO modifies the pdb structure and picks first chain'
	io = PDBIO()
	io.set_structure(structure[0]['A'])
	save_name = '%s_mod.pdb' %prot_name
	io.save('mod_pdb/'+save_name)

	'optional: show polypeptide chain of the modified pdb'
	if show_pp == 1:

		ppb = PPBuilder()
		for pp in ppb.build_peptides(structure[0]['A']):
			print(pp.get_sequence())
			print len(pp.get_sequence())
	else:
		pass

	if not check_residueNO == False:
		res_array_pos = check_residueNO - 1
		ppb = PPBuilder()
		for pp in ppb.build_peptides(structure[0]['A']):
			print pp.get_sequence()[res_array_pos]

	else:
		pass

def read_mumu(data_file):

	log_val = []
	'import data and split into a list by seperators'
	with open(data_file,'r') as f:
		
		for l in f:
			f_split_sep = l.split(':')
		
		'new list containing only log values'
		for i in range(len(f_split_sep)):
			log_val.append(f_split_sep[i].rstrip('\t\n').split(',')[:1][0])

		'remove header'	
		del log_val[0]
		
		'let us work in numpy arrays'
		log_array = np.asarray(log_val).astype(np.float)
		
		'saving to file'
		sv_name = "array_"+data_file
		np.savetxt(sv_name,log_array,delimiter=" ", fmt="%s")


print pdb_handling('2TRX_T77C_V91C',1,3)

#print read_mumu('2TRX_T77C_V91C_mumu.txt')















# def phaistos_handling(prot_name):
# 	current_path = os.getcwd()
# 	pdb_path = current_path+"/mod_pdb/"+prot_name+".pdb"


# 	os.chdir('/Users/thorn/phaistos/build2/bin')

# 	phaistos_command = './evaluate_observable --pdb-file '+pdb_path
# 	subprocess.call(phaistos_command)
	






