
def Mumu(phaistos_bin_dir,pdb_dir):
	import os
	import subprocess
	import numpy as np
	import pandas as pd
	import csv
	if phaistos_bin_dir == False:
		pass
	else:

		'phaistos need to excecuted from its directory'
		os.chdir(phaistos_bin_dir)

		'make a list of pdb files in the given folder'
		lsdir_nt = os.listdir(pdb_dir)
		'remove non-pdbs'
		lsdir = [x for x in lsdir_nt if ".pdb" in x is not True]

		
		log_p_sums = []
		log_p_id = []
		log_p_list = []

		'do the mumu for all files in the folder'
		for pp in range(len(lsdir)):

			cmd = "./evaluate_observable --pdb-file "+pdb_dir+lsdir[pp]+" --observable-mumu"

			direct_output = subprocess.check_output(cmd, shell=True)


			'we will clean it up a bit and put it into a list'
			log_val = []
			
			f_split_sep = direct_output.split(':')
			
			'new list containing only log values'
			for i in range(len(f_split_sep)):
				log_val.append(f_split_sep[i].rstrip('\t\n').split(',')[:1][0])

			'remove header'	
			del log_val[0]

				
			log_array = np.asarray(log_val).astype(np.float)

			'summing up log_values of residues'
			log_p_sums.append(sum(log_array))
			log_p_id.append(lsdir[pp])
			
			'save mumu-file for each protein'
			saving_dir = pdb_dir+"mumu/"
			if not os.path.exists(saving_dir):
				os.makedirs(saving_dir)
			saving_file = pdb_dir+"mumu/"+lsdir[pp].replace(".pdb", "")
			np.savetxt(saving_file,log_array,fmt="%s")
	
		delta_log_p = [-x + log_p_sums[0] for x in log_p_sums]


	
		'Quickly, let us add some free energy data - this simple task took me like 6 hours, silly.'
		'''ddg_file = "/Users/thorn/phaistos_proj/dd_G.txt"

		with open(ddg_file) as ddg:
			ddg_list = []
			for line in ddg:
				ddg_list.append(line.split("\r"))

		ddg_list_new = [0]
		for i in range(len(ddg_list[0])):
			ddg_list_new.append(ddg_list[0][i].replace("minus","-"))'''

		#free_energy = np.asarray(ddg_list_new).astype(np.float)
		#print len(free_energy)
		print len(log_p_id)
		print len(log_p_sums)
		print len(delta_log_p)
		np_id = np.array([log_p_id,log_p_sums,delta_log_p])

		#print np_id.T

		save_new = pdb_dir+"stability_data"
	
		np.savetxt(save_new,np_id.T,fmt="%s")


'Sys args magic'
import optparse
parser = optparse.OptionParser()

parser.add_option('-d', '--pha', help='Path to phaistos bin', dest='opt_phaistos_bin', default=False, action='store')
parser.add_option('-p', '--pdb', help='Path to pdb-file(s)', dest='opt_pdbs', default=False, action='store')

(opts, args) = parser.parse_args()


phaistos_bin_dir = opts.opt_phaistos_bin
pdb_dir = opts.opt_pdbs

Mumu(phaistos_bin_dir,pdb_dir)