
from Bio.PDB import *
import timeit

start = timeit.default_timer()

def print_full(x):
	import pandas as pd
	pd.set_option('display.max_rows', len(x))
	print(x)
	pd.reset_option('display.max_rows')

def read_protherm(input_csv):
	import pandas as pd
	from locale import *
	input_path = "df_input/"
	input_csv = input_path+input_csv
	#setlocale(LC_NUMERIC, 'de_DE.UTF-8')
	#For graphical purposes
	#pd.set_option('display.line_width', 5000) 
	#pd.set_option('display.max_columns', 60) 

	import numpy as np
	#import matplotlib.pyplot as plt

	pd.options.display.mpl_style = 'default'

	'Read the data'
	data = pd.read_csv(input_csv, sep=';', decimal = ',')

	'Clean up the data by defining what we want to remove'
	data = data[pd.notnull(data['ddG'])] #Remove rows with null data in ddG
	data = data.drop_duplicates('Mutation') #Remove rows with duplicates in Mutation
	
	'Let us handle commas like a boss'
	#data['ddG'] = data['ddG'].apply(atof)
	#data['ddG'] = data['ddG'].astype(float)

	data['Mutation'] = data['Mutation'].str.split(' ')
	
	'Return the dataframe'
	return data

def amino_index(PDB_name):
	'This function looks at a protein structure and retrieve its amino acids.'
	'Protein tructures often have amino acid chains that are truncated.'
	'This function correctly labels amino acids according to their position in the protein'
	from Bio.Seq import Seq
	from Bio import SeqIO
	from Bio.Alphabet import IUPAC
	from Bio.PDB import *
	io = PDBIO()
	import pandas as pd
	from Bio.PDB.DSSP import DSSP


	pdbl = PDBList()
	parser = PDBParser()
	ppb = PPBuilder()

	'If the structure hasnt been downloaed then it will - else parse it'
	structure = parser.get_structure(PDB_name,pdbl.retrieve_pdb_file(PDB_name))
	

	#io.set_structure(structure[0]['A'])

	'Choosing chain A'
	model = structure[0]
	chain = model['A']
	io.set_structure(chain)


	save_pdb = PDB_name+".pdb"
	io.save(save_pdb)
	
	'Get the sequence as a string'
	for pp in ppb.build_peptides(chain):
		seq = pp.get_sequence().lower()

	seq_list = []
	for i in range(len(seq)):
		seq_list.append(seq[i])

	'Return position of amino acids in the protein'
	pos_list = []
	for i in range(len(chain)):
		try:
			residue = chain[i]
			pos_list.append(i)

		except KeyError:
			pass

	'Dataframe with the amino acids position as its index'
	seq_pos = pd.DataFrame(seq_list, index=pos_list)

	return seq_pos

def save_mutated_chains(input_csv,PDB_name):
	import re
	import os

	def Save(PDB_name,mutation):
		true_i = mutation.index[i]
		
		'Save the chains with mutations as text files'
		save_dest = 'mutation_files/'+PDB_name+"/"
		if not os.path.exists(save_dest):
			os.makedirs(save_dest)

		mut_string = str(true_i)+"#"
		for j in range(len(mutation[true_i])):
			mut_string = mut_string+mutation[true_i][j]

		save_name = save_dest+mut_string
		text_file = open(save_name, "w")
		text_file.write(clean_chain.sum())
		text_file.close()

	mutation = read_protherm(input_csv)['Mutation']

	wt_amino_l = ['Non']
	mutation_pos_tag_l = ['Non']
	mutation_amino_tag_l = ['WILD']

	'Ok, we need to make a mutation but we also need to take care of double mutations!'
	for i in range(len(mutation)):
		clean_chain = amino_index(PDB_name)[0]

		true_i = mutation.index[i]
		'KeyError may occur - this is a work-around'
		try:
			if mutation[true_i][0] == "WILD":
				Save(PDB_name,mutation)

			else:
				'Divide by 3 to get number of mutations per chain'
				no_of_mut = len(mutation[true_i])/3

				'z is controls the index aka the position of each mutation if more than one are present'
				z = 1

				mutation_amino_tag = ""
				mutation_pos_tag = ""
				wt_amino_tag = ""
				for k in range(no_of_mut):

					mutation_pos = mutation[true_i][z]
					mutation_amino = mutation[true_i][z+1]
					wt_amino = mutation[true_i][z-1]

					'Remove some misplaced commas'
					mutation_pos = re.sub('[, ]', '', mutation_pos)
					mutation_amino = re.sub('[, ]', '', mutation_amino)

					mutation_amino_tag += mutation_amino
					wt_amino_tag += wt_amino
					mutation_pos_tag += mutation_pos+","

					

					clean_chain[int(mutation_pos)] = mutation_amino

					z =+ 4

					'Save the files'
					Save(PDB_name,mutation)

				mutation_pos_tag_l.append(mutation_pos_tag)
				mutation_amino_tag_l.append(mutation_amino_tag)
				wt_amino_l.append(wt_amino_tag)

		except KeyError:
				pass

	return [wt_amino_l, mutation_pos_tag_l, mutation_amino_tag_l]

def scwrl(PDB_name):
	import os 

	dir_name = "mutation_files/"+PDB_name
	lsdir = os.listdir(dir_name)

	'Remove hidden file'
	#lsdir.pop(0)

	dir_name_lvl1 = "mutation_files/"
	save_dest = dir_name_lvl1+"/"+PDB_name+"_mutated"+"/"
	if not os.path.exists(save_dest):
		os.makedirs(save_dest)

	for i in range(len(lsdir)):

		scw_path = "/Users/thorn/"
		pdb_path = PDB_name+".pdb"
		cmd = scw_path+"Scwrl4 -i "+pdb_path+" -o "+save_dest+"/"+lsdir[i]+"_mut"+".pdb -s " +dir_name+"/%s" %lsdir[i]
		os.system(cmd)

		print cmd

def mumu(phaistos_bin_dir,pdb_dir,input_csv):
	def add_dssp(pdb_file):
		'This is a small function that makes a dictionary with secondary structure information'

		from Bio.Seq import Seq
		from Bio import SeqIO
		from Bio.Alphabet import IUPAC
		
		io = PDBIO()
		import pandas as pd
		from Bio.PDB.DSSP import DSSP

		pdbl = PDBList()
		parser = PDBParser()
		ppb = PPBuilder()

		'If the structure hasnt been downloaed then it will - else parse it'

		structure = parser.get_structure(pdb_file,pdb_file)


		model = structure[0]
		chain = model['A']

		'path to the dssp excecutable:'
		dssp_exc = '/Users/thorn/dssp-2.2.1/mkdssp'
		dssp = DSSP(model, pdb_file, dssp_exc)

		sec_dict = {}
		for i in range(len(dssp)):
			a_key = dssp.keys()[i]
			index = a_key[1][1]
			sec_structure = dssp[a_key][1]
			sec_dict[index] = sec_structure

		return sec_dict




	import os
	import subprocess
	import numpy as np
	import pandas as pd
	import csv
	owd = os.getcwd()

	data_frame = read_protherm(input_csv)
	data_frame['log(p)'] = np.nan
	data_frame['d_log(p)'] = np.nan
	data_frame['dssp'] = np.nan

	if phaistos_bin_dir == False:
		pass
	else:

		'phaistos need to excecuted from its directory'
		os.chdir(phaistos_bin_dir)

		'make a list of pdb files in the given folder'
		lsdir_nt = os.listdir(pdb_dir)
		'remove non-pdbs'
		lsdir = [x for x in lsdir_nt if ".pdb" in x is not True]

		dssp_list = []
		'do the mumu for all files in the folder'
		for pp in range(len(lsdir)):

			pdb_path_ = pdb_dir+lsdir[pp]

			cmd = "./evaluate_observable --pdb-file "+pdb_path_+" --observable-mumu"

			direct_output = subprocess.check_output(cmd, shell=True)


			'we will clean it up a bit and put it into a list'
			log_val = []
			
			f_split_sep = direct_output.split(':')
			
			'new list containing only log values'
			for i in range(len(f_split_sep)):
				'converting -logp to logp'
				log_val.append(f_split_sep[i].rstrip('\t\n').split(',')[:1][0])

			'remove header'	
			del log_val[0]

			log_array = np.asarray(log_val).astype(np.float)

			'Converting -logp to logp'
			#log_array = np.multiply(log_array, -1.0)

			'Get the index of the mutation from the original dataframe by reading the number before #'
			index = lsdir[pp].partition("#")[0]
			sum_logp = sum(log_array)
			#print index

			pos = data_frame['Mutation'][int(index)]

			if len(pos) == 1:
				dssp_info = None
			else:
				a_index = int(pos[1])
				dssp_info = add_dssp(pdb_path_)[a_index]

			dssp_list.append(str(dssp_info))

			data_frame['log(p)'][int(index)] = sum_logp
			data_frame['d_log(p)'][int(index)] = data_frame['log(p)'][0]-sum_logp
	
	data_frame['dssp'] = dssp_list


	return data_frame





import pandas as pd
import os

#print read_protherm('1BPI.csv')

#1YCC 451C' '1LZ1'
protein_list = ['1FTG','1CSP','1VQB','1BNI','1ROP','1BPI','2RN2','1STN','1BVC','2LZM','1IGV','3SSI','1CYO','2TRX','1PGA','1C9O']
#protein_list = ['1FTG']

PDB_name_list = protein_list
input_csv_list = []
for i in range(len(protein_list)):
	input_csv_list.append(protein_list[i]+".csv")

owd = os.getcwd()

for i in range(len(input_csv_list)):
	os.chdir(owd)
	input_csv = input_csv_list[i]
	PDB_name = PDB_name_list[i]


	phaistos_bin_dir = "/Users/thorn/phaistos/build/bin/"
	pdb_dir = "/Users/thorn/phaistos_proj/project2/mutation_files/"+PDB_name+"_mutated/"
	


	#print print_full(read_protherm(input_csv)['Mutation'])
	df_mut = save_mutated_chains(input_csv,PDB_name)
	scwrl(PDB_name)
	data_frame = mumu(phaistos_bin_dir,pdb_dir,input_csv)

	data_frame['WT_AA'] = df_mut[0]
	data_frame['Position'] = df_mut[1]
	data_frame['Mut_AA'] = df_mut[2]


	

	save_path = "/Users/thorn/phaistos_proj/project2/df_output/"
	save_name = PDB_name+"_data_frame.csv"

	if not os.path.exists(save_path):
		os.makedirs(save_path)
	data_frame.to_csv(save_path+save_name, sep=';')

	print print_full(data_frame)

stop = timeit.default_timer()

print stop - start




