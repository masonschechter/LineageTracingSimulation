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
(["A", "C", "G", "T", "T"],["A", "C", "G", "T", "T", "C", "G", "A"])
]

illumina_adapter = list('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')

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

def barcode(well, illumina_adapter, fwd_barcode, rev_barcode):
	## Generates complete sequence to include illumina adapter and each barcode ##
	fwd_seqs = []
	rev_seqs = []
	for cell in well:
		umi = list(np.random.choice(["G", "C", "A", "T"], size=20))
		fwd_seq = illumina_adapter + fwd_barcode + cell + umi + reverse_complement(rev_barcode)
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



well_1 = []
well_1_1 = []
well_1_2 = []
well_1_1_1 = []
well_1_1_2 = []
well_1_2_1 = []
well_1_2_2 = []

fwd_sequences = []
rev_sequences = []

for x in range(0, initial_cell_count):
	well_1.append(list(WT_seq))

well_1 = cell_cycle(well_1, cut_freq, insert_freq, del_freq)
well_1_1, well_1_2 = split_cells(well_1)

well_1_1 = cell_cycle(well_1_1, cut_freq, insert_freq, del_freq)
well_1_2 = cell_cycle(well_1_2, cut_freq, insert_freq, del_freq)
well_1_1_1, well_1_1_2 = split_cells(well_1_1)
well_1_2_1, well_1_2_2 = split_cells(well_1_2)

well_1_1_1 = cell_cycle(well_1_1_1, cut_freq, insert_freq, del_freq)
well_1_1_2 = cell_cycle(well_1_1_2, cut_freq, insert_freq, del_freq)
well_1_2_1 = cell_cycle(well_1_2_1, cut_freq, insert_freq, del_freq)
well_1_2_2 = cell_cycle(well_1_2_2, cut_freq, insert_freq, del_freq)


well_1_1_1 = seq_filter(well_1_1_1)
well_1_1_2 = seq_filter(well_1_1_2)
well_1_2_1 = seq_filter(well_1_2_1)
well_1_2_2 = seq_filter(well_1_2_2)

well_1_1_1_fwd, well_1_1_1_rev = barcode(well_1_1_1, illumina_adapter, fwd_rev_barcodes[0][0], fwd_rev_barcodes[0][1])
well_1_1_2_fwd, well_1_1_2_rev = barcode(well_1_1_2, illumina_adapter, fwd_rev_barcodes[1][0], fwd_rev_barcodes[1][1])
well_1_2_1_fwd, well_1_2_1_rev = barcode(well_1_2_1, illumina_adapter, fwd_rev_barcodes[2][0], fwd_rev_barcodes[2][1])
well_1_2_2_fwd, well_1_2_2_rev = barcode(well_1_2_2, illumina_adapter, fwd_rev_barcodes[3][0], fwd_rev_barcodes[3][1])

well_1_1_1_fwd = amplify(well_1_1_1_fwd)
well_1_1_2_fwd = amplify(well_1_1_2_fwd)
well_1_2_1_fwd = amplify(well_1_2_1_fwd)
well_1_2_2_fwd = amplify(well_1_2_2_fwd)

well_1_1_1_rev = amplify(well_1_1_1_rev)
well_1_1_2_rev = amplify(well_1_1_2_rev)
well_1_2_1_rev = amplify(well_1_2_1_rev)
well_1_2_2_rev = amplify(well_1_2_2_rev)

fwd_sequences = convert_to_str(well_1_1_1_fwd + well_1_1_2_fwd + well_1_2_1_fwd + well_1_2_2_fwd)
rev_sequences = convert_to_str(well_1_1_1_rev + well_1_1_2_rev + well_1_2_1_rev + well_1_2_2_rev)

