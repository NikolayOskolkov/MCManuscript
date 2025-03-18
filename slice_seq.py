# This script splits a fasta-file into 60bp reads in a jumping window approach. 
# Run this script as: python slice_seq.py
import os

os.chdir('/home/nikolay/WABI/T_van_der_Valk/sliced_seq/')

n_seqs = 0
with open('output_sliced_seqs/GTDB_sliced_seqs.fna', 'w') as fout:
	
	fna_files = os.listdir('fasta')
	for i in fna_files:
		print('Working with file ' + str(i))
		
		with open('fasta/' + str(i),'r') as fin:
			for line in fin:
				if len(line)==61 and 'x' and 'N' not in line:
					n_seqs = n_seqs + 1
					if(n_seqs % 10000 == 0):
						print(n_seqs)
					fout.write('>' + str(i) + '_' + str(n_seqs) + '\n')
					fout.write(line)
