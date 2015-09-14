#### Author ####
# Theis Hjalte Thorn Jakobsen
# Please cite if you decide to use this script for your data
# Questions may be directed to thornjakobsen@gmail.com

#### Imports & dependicies ####
import pandas as pd

class sortingClass(object):

	def __init__(self, csvfile):
		self.csvfile = csvfile

	# This will print a complete dataframe (normally pandas dataframes are shortened).
	def print_full(self):
		pd.set_option('display.max_rows', len(self.df))
		return(self.df)
		pd.reset_option('display.max_rows')

	def organize_data(self):
		# This function reads a csv file with protherm data.
		# Null values in target columns are removed.
		# Data is grouped by the proteins wild-type PDB name.
		self.df = pd.read_csv(self.csvfile, sep=';', decimal = '.').sort()
		self.df = self.df.dropna(subset=(['ddG','PDB_wild','Mutation']))

		# Commas are removed and replaced by a dot.
		no_comma = lambda x: float(x.replace(',','.'))
		self.df['ddG'] = self.df['ddG'].apply(no_comma)

		# The is all the proteins we got that have ddG values.
		unique_pdb = self.df['PDB_wild'].unique()

		# The data is grouped by a column name.
		self.df_grp = self.df.groupby(['PDB_wild'])

		# Remove certain duplicates within subsets.
		self.df_dict = {}
		for i, k in self.df_grp:
			k.drop_duplicates(subset=['Mutation'])
			self.df_dict[i] = k

		return self.df_dict

	# Counter for counting amount of proteins and mutations left
	def counter(self):
		self.sum_proteins = 0
		self.sum_mutations = 0
		for i in self.df_dict.keys():
			self.sum_proteins += 1
			self.sum_mutations += len(self.df_dict[i])

		return self.sum_proteins, self.sum_mutations
