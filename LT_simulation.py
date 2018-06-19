import numpy as np

#WT_seq = 'CTTGTGGAAAGGACGAAACACCGGCCCAGACTGAGCACGTGAGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTAAGCTTGGGCCGCTCGAGGTACCTCTCTACATA'
# fwd_read_file = '/data/users/mscheche/simulation/fwd_reads.fastq'
# rev_read_file = '/data/users/mscheche/simulation/rev_reads.fastq'
#WT seq from March 29 run

WT_seq = 'ATGCGATCAGATGCTACGTGCATGACAGATCA'
fwd_read_file = 'C:/Users/Mason/Documents/simulation_output/fwd_reads.fastq'
rev_read_file = 'C:/Users/Mason/Documents/simulation_output/rev_reads.fastq'

fwd_rev_barcodes = [
(["A", "G", "A", "A", "G"],["A", "G", "A", "A", "G", "T", "C", "T", "A", "A"]),
(["T", "A", "C", "C", "T"],["T", "A", "C", "C", "T", "T", "G", "C", "T", "G", "G"]), 
(["G", "T", "A", "G", "G"],["G", "T", "A", "G", "G", "C", "A", "A", "T"]),
(["A", "C", "G", "T", "T"],["A", "C", "G", "T", "T", "C", "G", "A"]),
(["C", "G", "T", "A", "A"],["G", "A", "A", "A", "T", "C", "C"])
(["G", "G", "A", "C", "T"],["C", "A", "G", "C", "G", "T", "A", "G"])
(["T", "A", "C", "C", "T"],["A", "G", "A", "A", "G", "T", "C", "T", "A", "A"])
(["G", "T", "A", "G", "G"],["T", "A", "C", "C", "T", "T", "G", "C", "T", "G", "G"])
(["A", "C", "G", "T", "T"],["G", "T", "A", "G", "G", "C", "A", "A", "T"]),
(["C", "G", "T", "A", "A"],["A", "C", "G", "T", "T", "C", "G", "A"]),
(["G", "G", "A", "C", "T"],["G", "A", "A", "A", "T", "C", "C"]),
(["G", "T", "A", "G", "G"],["A", "G", "A", "A", "G", "T", "C", "T", "A", "A"]),
(["A", "C", "G", "T", "T"],["T", "A", "C", "C", "T", "T", "G", "C", "T", "G", "G"]),
(["C", "G", "T", "A", "A"],["G", "T", "A", "G", "G", "C", "A", "A", "T"]),
(["G", "G", "A", "C", "T"],["A", "C", "G", "T", "T", "C", "G", "A"]),
(["A", "G", "A", "A", "G"],["C", "A", "G", "C", "G", "T", "A", "G"])
]

illumina_adapter = [list('ACACTCTTTCCCTACACGACGCTCTTCCGATCT'), list('TTCAGACGTGTGCTCTTCCGATCT')]

quality_scores = ["<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", ":", ";"]

##############Parameters###############
initial_cell_count = 100
num_of_generations = 3
cut_freq = .40 
insert_freq = .80
del_freq = .20
insertion_len_range = [0, 1, 2, 3, 4, 5]
insertion_len_freq = [.15, .20, .20, .30, .10, .05]
nucleotide_bias = [.35, .30, .20, .15]
insert_pos = -8 ##39/40
########################################

def does_it_cut(cut_freq):
	## Binary output to determine if a cut occurs, based on observed cutting frequency ##
	cut = np.random.choice((1,0), p=[cut_freq, 1-cut_freq])
	return cut

def insertion_or_deletion(insert_freq, del_freq):
	## Binary output if edit is an insertion or deletion, based on observed insertion and deletion frequencies ##
	indel = np.random.choice((1,0), p=[insert_freq, del_freq])
	return indel

def split_cells(well_to_split):
	## Simulates splitting of cells, roughly but not precisely 50-50 ##
	next_well_1 = []
	next_well_2 = []
	for cell in well_to_split:
		split = np.random.choice((1,0))
		if split:
			next_well_1.append(cell)
		else:
			next_well_2.append(cell)
	return next_well_1, next_well_2

def which_nucleotide():
	## Returns "random" nucleotide based on Tdt nucleotide bias ##
	nucleotides = ["G", "C", "A", "T"]
	inserted_nucleotide = np.random.choice(nucleotides, p=nucleotide_bias)
	return inserted_nucleotide

def cell_division(old_well):
	## Simulates cell division ##
	new_well = []
	for x in range(0, len(old_well)):
		new_well.extend([old_well[x], old_well[x]])
	return new_well

def cut_indel(well, cut_freq, insert_freq, del_freq):
	## Simulates an editing event based on frequency of cutting, insertions, and deletions ##
	new_well = []
	for x in range(0, len(well)):
		cell = well[x]
		cut = does_it_cut(cut_freq)
		if not cut:
			new_well.append(cell)
			continue
		insertion = insertion_or_deletion(insert_freq, del_freq)
		if insertion:
			insertion_len = np.random.choice(insertion_len_range, p=insertion_len_freq)
			for x in range(0, insertion_len):
				nucleotide = which_nucleotide()
				cell.insert(insert_pos, nucleotide)
			new_well.append(cell)
		else:
			for x in range(2, -1, -1):
				del cell[insert_pos-x]
			new_well.append(cell)
	return new_well

def reverse_complement(seq):
	## Generates reverse complement of target sequence ##
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	reverse_seq = []
	for base in seq[::-1]:
		reverse_seq.append(complement[base])
	return reverse_seq

def barcode(well, fwd_illumina_adapter, rev_illumina_adapter, fwd_barcode, rev_barcode):
	## Generates complete sequence to include illumina adapter and each barcode ##
	fwd_seqs = []
	rev_seqs = []
	for cell in well:
		umi = list(np.random.choice(["G", "C", "A", "T"], size=20))
		fwd_seq = fwd_illumina_adapter + fwd_barcode + cell + umi + reverse_complement(rev_barcode) + reverse_complement(rev_illumina_adapter)
		rev_seq = reverse_complement(fwd_seq)
		fwd_seqs.append(fwd_seq)
		rev_seqs.append(rev_seq)
	return fwd_seqs, rev_seqs

def cell_cycle(well, cut_freq, insert_freq, del_freq):
	## Ideal cell cycle with a possible editing event before cell division ##
	well = cell_division(cut_indel(well, cut_freq, insert_freq, del_freq))
	return well

def generate_fastq(read_file, seq_list):
	## Writes out to a text file in fastq format ##
	with open(read_file, 'w') as file:
		date = datetime.now().strftime('%H:%M:%S')
		for seq in seqs:
			file.write('@:' + date + '\n')
			file.write(seq + '\n')
			file.write('+\n')
			for x in range(0, len(seq)):
				file.write(str(np.random.choice(quality_scores)))
			file.write('\n')

def pcr(seq):
	## Simulates NGS error rate during amplification ##
	amplicon = []
	for x in range(0, len(seq)):
		mutation = np.random.choice([1,0], p=[.001, .999]) ##NGS mutation rate
		if not mutation:
			amplicon.append(seq[x])
		else:
			amplicon.append(np.random.choice(["G", "C", "A", "T"]))
	return amplicon

def amplify(seq_list):
	## Simulates NGS amplification over a range ##
	amplified_seqs = []
	for seq in seq_list:
		num = np.random.choice(range(20, 40))
		for x in range(0, num):
			amplified_seqs.append(pcr(seq))
	return amplified_seqs


def convert_to_str(seq_list):
	## Seqs as list of characters => seqs as strings to write out to file ##
	seqs_as_str = []
	for seq in seq_list:
		seqs_as_str.append("".join(seq))
	return seqs_as_str

def seq_filter(seq_list):
	## Number of unique molecules collected in sample ##
	chosen_seqs = []
	for seq in seq_list:
		keep = np.random.choice([1,0], p=[.2, .8])
		if keep:
			chosen_seqs.append(seq)
	return chosen_seqs


wells = []
for num in range(0, (2**num_of_generations)-1):
	wells.append([])

for x in range(0, initial_cell_count):
	wells[0].append(list(WT_seq))

for i, well in enumerate(wells):
	if i == 0:
		wells[0] = cell_cycle(wells[0], cut_freq, insert_freq, del_freq)
		wells[1], wells[2] = split_cells(well)
	else:
		wells[i] = cell_cycle(well, cut_freq, insert_freq, del_freq)
		if ((2*i)+2) <= len(wells):
			wells[(2*i)+1], wells[(2*i)+2] = split_cells(wells[i])


last_generation = []
for i in range(0, (num_of_generations**2)//2):
	last_generation.append([])

for i in range(int(-1*(np.rint(len(wells)/2))), 0):
	last_generation[i] = seq_filter(wells[i])

fwd_barcoded_seqs = []
rev_barcoded_seqs = []
for i, (x, y) in enumerate(fwd_rev_barcodes):
	if i < len(last_generation):
		print(i, (x,y))
		fwd_well, rev_well = barcode(last_generation[i], illumina_adapter[0], illumina_adapter[1], x, y)
		fwd_barcoded_seqs += amplify(fwd_well)
		rev_barcoded_seqs += amplify(rev_well)

fwd_reads = convert_to_str(fwd_barcoded_seqs)
rev_reads = convert_to_str(rev_barcoded_seqs)

generate_fastq(fwd_read_file, fwd_reads)
generate_fastq(rev_read_file, rev_reads)