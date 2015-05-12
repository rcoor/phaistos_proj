#### Author ####
# Theis Hjalte Thorn Jakobsen
# Please cite if you decide to use this script for your data
# Questions may be directed to thornjakobsen@gmail.com

#### Imports & dependicies ####
import pandas as pd


# This will print the complete data-frame.
def print_full(x):
	import pandas as pd
	pd.set_option('display.max_rows', len(x))
	print(x)
	pd.reset_option('display.max_rows')

def protherm_read_and_group(protherm_csv):
	
	# This function reads a csv file with protherm data.
	# Null values in target columns are removed.
	# Data is grouped by the proteins wild-type PDB name.
	input_csv = protherm_csv

	df = pd.read_csv(input_csv, sep=';', decimal = '.').sort()
	
	df = df.dropna(subset=(['ddG','PDB_wild','Mutation']))
	
def no_comma(x):
		x = float(x.replace(',','.'))
		return x

	# ddG data is converted to float.
	df[['ddG','T']] = df[['ddG','T']].apply(no_comma)

	# The is all the proteins we got that have ddG values.
	unique_pdb = df['PDB_wild'].unique()

	# The data is grouped by a column name.
	df_grp = df.groupby(['PDB_wild'])
	
	# Remove certain duplicates within subsets.
	df_dict = {}
	for i, k in df_grp:
		
		k.drop_duplicates(subset=['Mutation'])

		df_dict[i] = k

	return df_dict
		

print protherm_read_and_group('alldata.csv')


