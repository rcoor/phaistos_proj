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

PDB_name = '1STN.pdb'
structure = parser.get_structure(PDB_name,PDB_name)


model = structure[0]
chain = model['A']

dssp_exc = '/Users/thorn/dssp-2.2.1/mkdssp'
dssp = DSSP(model, '1STN.pdb', dssp_exc)

#print dssp['res_id']

sec_dict = {}
for i in range(len(dssp)):
	a_key = dssp.keys()[i]
	index = a_key[1][1]
	sec_structure = dssp[a_key][1]

	sec_dict[index] = sec_structure

print sec_dict[10]



	

