# Here we import needed scripts.
from sorting import sortingClass

from get_protein import getPDB_mutate

import os
# from mutate import save_mutated_chains



class mutate(object):

    def __init__(self, datafile):

        self.datafile = datafile
        # Using the sortingClass for cleaning and grouping protein data.
        # The organize_data creates a dictionary with proteins names.
        self.readdata = sortingClass(self.datafile)
        self.proteinsData = self.readdata.organize_data()
    def strings(self):

        # Calls organize_data to initialize protein dictionary.

        self.readdata.counter()

        # amino_index gets a list of the placement of all amino acids in the chain of a given protein.
        # for i in proteinsData.keys():
        #     key = proteinsData.keys()[i]
        #     PDBClass = getPDB_mutate(key, proteinsData[key])

        for i in range(len(self.proteinsData.keys())):
            try:

                key = self.proteinsData.keys()[i]
                #print proteinsData[key]
                PDBClass = getPDB_mutate(key, self.proteinsData[key])

                PDBClass.amino_index()
                PDBClass.make_mutation()
            except (KeyError, ValueError):
                pass

# Let's start mutating strings.
mutate('alldata.csv').strings()
# Counting amount of mutations.
cpt = sum([len(files) for r, d, files in os.walk("mutation_files")])
print cpt
