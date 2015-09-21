# Here we import needed scripts.
from sorting import sortingClass

from get_protein import getPDB_mutate
# from mutate import save_mutated_chains

# Using the sortingClass for cleaning and grouping protein data.
# The organize_data creates a dictionary with proteins names.
readdata = sortingClass('alldata.csv')

# Calls organize_data to initialize protein dictionary.
proteinsData = readdata.organize_data()
readdata.counter()

# amino_index gets a list of the placement of all amino acids in the chain of a given protein.
# for i in proteinsData.keys():
#     key = proteinsData.keys()[i]
#     PDBClass = getPDB_mutate(key, proteinsData[key])

i = 5
key = proteinsData.keys()[i]
PDBClass = getPDB_mutate(key, proteinsData[key])

PDBClass.amino_index()
print PDBClass.make_mutation()
