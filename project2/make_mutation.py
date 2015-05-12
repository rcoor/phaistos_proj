def make_mutation():
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	import read_protherm

	return read_protherm('Barnase.csv')





print make_mutation()