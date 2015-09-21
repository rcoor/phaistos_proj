from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.PDB import *
io = PDBIO()
import pandas as pd
from Bio.PDB.DSSP import DSSP
import os
import re
class getPDB_mutate(object):

    def __init__(self, PDB_name, df):
        self.PDB_name = PDB_name
        self.df = df
        'This function looks at a protein structure and retrieve its amino acids.'
        'Protein tructures often have amino acid chains that are truncated.'
        'This function correctly labels amino acids according to their position in the protein'
        self.pdbl= PDBList()
        self.parser= PDBParser()
        self.ppb= PPBuilder()

        'If the structure hasnt been downloaded then it will - else parse it'
        structure = self.parser.get_structure(self.PDB_name,self.pdbl.retrieve_pdb_file(self.PDB_name))

        #io.set_structure(structure[0]['A'])
        'Choosing chain A'
        model = structure[0]
        self.chain = model['A']
        io.set_structure(self.chain)

        if not os.path.exists("wildtypes"):
            os.makedirs("wildtypes")

        save_pdb = "wildtypes/"+self.PDB_name+".pdb"
        io.save(save_pdb)

    def amino_index(self):

        'Get the sequence as a string'
        for pp in self.ppb.build_peptides(self.chain):
            seq = pp.get_sequence().lower()

        seq_list = []
        for i in range(len(seq)):
            seq_list.append(seq[i])

        'Return position of amino acids in the protein'
        pos_list = []
        for i in range(len(self.chain)):
            try:
                residue = self.chain[i]
                pos_list.append(i)

            except KeyError:
                pass

        'Dataframe with the amino acids position as its index'
        self.sequencePositions = pd.DataFrame(seq_list, index=pos_list)[0]
        return self.sequencePositions

    def make_mutation(self):
        self.sequencePositions

        for i in range(len(self.df['Mutation'])):
            print i
            true_i = self.df['Mutation'].index[i]
            print true_i

            try:
                if self.df['Mutation'][0] == "WILD":
                    print self.df['Mutation'][true_i][0]
                    save_dest = 'mutation_files/'+self.PDB_name+'/'
                    if not os.path.exists(save_dest):
                        os.makedirs(save_dest)
                else:
                    pass
            except KeyError:
                pass
    #             if self.df['Mutation'][0] == "WILD":
    #                 print self.mutation[true_i][0]
    #                 'Save the chains with mutations as text files'
    #                 save_dest = 'mutation_files/'+self.PDB_name+"/"
    #                 if not os.path.exists(save_dest):
    #                     os.makedirs(save_dest)
    #
    #                 mut_string = str(true_i)+"#"
    #                 for j in range(len(self.mutation[true_i])):
    #                     mut_string = mut_string+self.mutation[true_i][j]
    #
    #                 save_name = save_dest+mut_string
    #                 text_file = open(save_name, "w")
    #                 text_file.write(self.sequencePositions.sum())
    #                 text_file.close()
    #
    #             else:
    #                 'Divide by 3 to get number of mutations per chain'
    #                 no_of_mut = len(mutation[true_i])/3
    #
    #         except KeyError:
    #           pass
