from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import pickle as pkl
import argparse
import sys
sys.path.append('library')
from library.assemble_sentences import *
import time


#############################################################
########################  Parameters  #######################
#############################################################
parser = argparse.ArgumentParser(description="HOTSPOT is a learning-based tool designed for plasmid host prediction. Its backbone is a phylogenetic tree of plasmids' hosts (bacteria) from phylum to species. The top-down tree search can accurately predict the hosts' taxonomic labels by incorporating the state-of-the-art language model, Transformer, in each node’s taxon classifier. To use HOTSPOT, you only need to input plasmid sequences (complete or segmented) into the program.")
parser.add_argument('--contigs', help='FASTA file of the input sequences (one or more contigs in a single FASTA file, default test_contigs.fa)',  default = 'test_contigs.fa')
parser.add_argument('--len', help='minimum length (bp) of contigs for length filtering (default 1500)', type=int, default=1500)
parser.add_argument('--threads', help='number of threads to use (default 8)', type=int, default=8)
parser.add_argument('--dbdir', help='database directory (default database/)',  default = 'database/')
parser.add_argument('--midfolder', help='folder to store the intermediate files from preprocessing (default temporary_files/)', type=str, default='temporary_files/')
inputs = parser.parse_args()


#############################################################
######################  Check folders  ######################
#############################################################
def help_info():
    print('')
    print("""The usage of preprocessing.py:
            [-h, --help]          Show the help message and exit
            [--contigs INPUT_FA]  FASTA file of the input sequences (one or more contigs in a single FASTA file, default test_contigs.fa)
            [--len MINIMUM_LEN]   Minimum length (bp) of contigs for length filtering (default 1500)
            [--threads NUM]       Number of threads to use (default 8)
            [--dbdir DR]          Database directory (default database/)
            [--midfolder DIR]     Folder to store the intermediate files from preprocessing (default temporary_files/)""")


out_fn = inputs.midfolder

if not os.path.isdir(out_fn):
    os.makedirs(out_fn)

db_dir = inputs.dbdir
if not os.path.exists(db_dir):
    print(f"""Database directory "{db_dir}" unavailable. Please use the option "--dbdir" to specify the database path or place the database folder under HOTSPOT's main directory ("HOTSPOT/database/").""")
    help_info()
    exit(1)

if(os.path.exists(inputs.contigs)==False):
    print(f'Input FASTA file not found. Please use the option "--contigs" to specify the input plasmid contigs')
    help_info()
    exit(1)


#############################################################
##################  Filter short contigs  ###################
#############################################################
rec = []
length_dict = {}
for record in SeqIO.parse(inputs.contigs, 'fasta'):
    length_dict[record.id] = len(record.seq)
    if len(record.seq) >= inputs.len:
        rec.append(record)
if(len(rec)==0):
    print(f'No contig left after length filtering!')
    exit(1)
SeqIO.write(rec, f'{out_fn}/filtered_contigs.fa', 'fasta')
#pkl.dump(length_dict, open(f'{out_fn}/length.dict', 'wb'))


#############################################################
#####################  Preprocessing  #######################
#############################################################
start = time.time()
threads = inputs.threads
preprocessing(f'{out_fn}/filtered_contigs.fa', db_dir, out_fn, threads)
sentence = assemble_sentences(out_fn, db_dir)
if(len(sentence)==0):
    print(f'The encoded sentence matrix are empty! (no sufficient feature for host prediction)')
    exit(1)


end = time.time()
print(f'The total time for preprocessing.py is {end-start}s.')
