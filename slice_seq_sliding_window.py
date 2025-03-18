# This script splits a fasta-file into 60bp reads in a sliding window manner with 10bp step. 
# Run this script as: /home/nikolay/miniconda3/bin/python slice_seq_sliding_window.py

import os
from Bio import SeqIO

os.chdir('/home/nikolay/WABI/T_van_der_Valk/sliced_seq_sliding_window/')

def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0, seqlen, step):
        j = seqlen if i + win > seqlen else i + win
        yield seq[i:j]
        if j == seqlen: break

n_seqs = 0
window_size = 60
window_step = 10
with open('output_sliced_seqs/GTDB_sliced_seqs_sliding_window.fna', 'w') as fout:
	
	fna_files = os.listdir('fasta')	
	for i in fna_files:
		print('Working with file ' + str(i))
		
		fasta_sequences = SeqIO.parse('fasta/' + str(i), 'fasta')
		for fasta in fasta_sequences:
			
			name, sequence = fasta.id, str(fasta.seq)
			for subseq in chunks(sequence, window_size, window_step):
				if len(subseq) == window_size and 'x' not in str(subseq):
					n_seqs = n_seqs + 1
					if(n_seqs % 10000 == 0):
						print(n_seqs)
					fout.write('>' + str(i) + '_' + str(n_seqs) + '\n')
					fout.write(subseq + '\n')
		
