# Here we import needed scripts.
from sorting import sortingClass
from get_protein import getPDB

# Using the sortingClass for cleaning and grouping protein data.
# The organize_data creates a dictionary with proteins names.
readdata = sortingClass('alldata.csv')

# Calls organize_data to initialize protein dictionary.
proteinsData = readdata.organize_data()
print readdata.counter()

# amino_index gets a list of the placement of all amino acids in the chain of a given protein.
print getPDB(proteinsData.keys()[0]).amino_index()
