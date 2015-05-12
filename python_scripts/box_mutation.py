


def make_mutation(csv_file,save_folder,tru):
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	import re
	import csv
	import os
	
	'Import CSV list of mutations with wt sequence as header'
	mutations_list = csv.reader(open(csv_file,'rU'), delimiter=";")
	mutations = []
	for row in mutations_list:
		mutations.append(row[0])

	'split into sequence and mutations'
	prot_seq = mutations[0].lower()
	mutations.pop(0)

	'Creates a folder for data-output'
	directory = save_folder
	if not os.path.exists(directory):
		os.makedirs(directory)

	'Makes a table with all mutations with their matching sequence'
	table_of_mutations = [("aaawildtype",str(Seq(prot_seq)[int(tru):]))]
	i = 0
	for i in range(len(mutations)):
		
		my_seq_un = Seq(prot_seq)
		my_seq = my_seq_un.tomutable()

		if mutations == 'false':
			pass
		else:
			match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", mutations[i], re.I)
			if match:
				position_id = match.groups()
	    		
	    		position = position_id[1]

	    		my_seq[int(position)-1] = position_id[2]
	    		
	    	
	    		if mutations[i]:
	    			table_of_mutations.append((mutations[i],my_seq[int(tru):]))

	    			'save as single files, used for Scwrl'
	    			with open(os.path.join(directory, str(mutations[i])), 'wb') as temp_file:
	    				temp_file.write(str(my_seq[int(tru):]))

	with open(os.path.join(directory, "aaawildtype"), 'wb') as temp_file:
	    temp_file.write(str(Seq(prot_seq)[int(tru):]))


	'Saving the list as a file'
	import csv
	with open(os.path.join(directory,'sum_all_mutations.csv'),'w') as out:
		csv_out = csv.writer(out)
		csv_out.writerow(['Mutation','Sequence'])
		for row in table_of_mutations:
			csv_out.writerow(row)


def Scwrl(scw,dir_name,pdb_path):
	

	if scw == False:
		pass
	else:
		if pdb_path == False:
			pass
		else:
			import os

			scw_path = scw

			lsdir = os.listdir(dir_name)

			lsdir.remove("sum_all_mutations.csv")

			save_dest = dir_name+"pdbs_mut"
			if not os.path.exists(save_dest):
				os.makedirs(save_dest)

			for j in range(len(lsdir)):

				print lsdir[j]
				cmd = scw_path+"Scwrl4 -i "+pdb_path+" -o "+save_dest+"/MUT_"+lsdir[j]+".pdb -s" +dir_name+"/%s" %lsdir[j]
				os.system(cmd)

				print cmd

		





import subprocess as sub
import optparse

'Sys args magic'	
parser = optparse.OptionParser()

parser.add_option('-c', '--csv', help='Specify a csv file', dest='opt_csv', default=False, action='store')
parser.add_option('-d', '--dir', help='Specify a saving directory', dest='opt_dir', default=False, action='store')
parser.add_option('-s', '--scw', help='Mutate with Scwrl, specify path to Scwrl', dest='opt_scw', default=False, action='store')
parser.add_option('-t', '--tru', help='N-terminal truncation, default = 0', dest='opt_tru', default=0, action='store')
parser.add_option('-p', '--pdb', help='Input pdb file', dest='opt_pdb', default=False, action='store')
(opts, args) = parser.parse_args()

if opts.opt_csv == False:
	print "specify an input file with -c <file>"
else:
		csv_file = opts.opt_csv

if opts.opt_dir == False:
	print "specify saving directory -d <file>"
else:
		dir_name = opts.opt_dir

if opts.opt_scw == False:
	scw = False
else:
	scw = opts.opt_scw

tru = opts.opt_tru

pdb_path = opts.opt_pdb

'''Gotta run 'em all'''
make_mutation(csv_file,dir_name,tru)
Scwrl(scw,dir_name,pdb_path)

