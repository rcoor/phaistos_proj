
import os
class mutatePDB(object):

	def __init__(self, mut_dir):
		self.mut_dir = mut_dir


	def get_pdb_names(self):
		self.mut_dir
		self.dir_content = os.listdir(self.mut_dir)
		self.dir_content = filter(lambda a: len(a) == 4, self.dir_content)

	def scwrl(self):

		for self.PDB_name in self.dir_content:

			dir_name = self.mut_dir+"/"+self.PDB_name
			lsdir = os.listdir(dir_name)

			'Remove hidden file'
			#lsdir.pop(0)

			dir_name_lvl1 = "mutation_files/"
			save_dest = dir_name_lvl1+"/"+self.PDB_name+"_mutated"+"/"
			if not os.path.exists(save_dest):
				os.makedirs(save_dest)

			for i in range(len(lsdir)):

				scw_path = "/Users/thorn/"
				pdb_path = "wildtypes/"+self.PDB_name+".pdb"
				cmd = scw_path+"Scwrl4 -i "+pdb_path+" -o "+save_dest+"/"+lsdir[i]+"_mut"+".pdb -s " +dir_name+"/%s" %lsdir[i]
				os.system(cmd)

				print cmd
