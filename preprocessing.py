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
parser = argparse.ArgumentParser(description="""HOTSPOT is a learning-based tool for plasmid host prediction. Its backbone is a phylogenetic tree of the plasmid hosts (bacteria) from phylum to species. By incorporating the state-of-the-art language model, Transformer, in each node’s taxon classifier, the top-down tree search can accurately predict the host taxonomy for the input plasmid contigs. There are totally 115 taxon classifiers, each corresponding to a node with more than one child node. To use HOTSPOT, you only need to input complete plasmids or plasmid contigs assembled from metagenomic data into the program. 
                                                """)
parser.add_argument('--contigs', help='FASTA file of the inputs (one or more contigs, default test_contigs.fa)',  default = 'test_contigs.fa')
parser.add_argument('--len', help='minimum length (bp) of contigs for length filter (default 1500)', type=int, default=1500)
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
            [--contigs INPUT_FA]  FASTA file of the inputs (one or more contigs, default test_contigs.fa)
            [--len MINIMUM_LEN]   Minimum length (bp) of contigs for length filter (default 1500)
            [--threads NUM]       Number of threads to use (default 8)
            [--dbdir DR]          Database directory (default database/)
            [--midfolder DIR]     Folder to store the intermediate files from preprocessing (default temporary_files/)
    """)


out_fn = inputs.midfolder

if not os.path.isdir(out_fn):
    os.makedirs(out_fn)

db_dir = inputs.dbdir
if not os.path.exists(db_dir):
    print(f'Database directory "{db_dir}" missing or unreadable. Please use option "--dbdir" to specify the database path or place the database folder under the main directory of HOTSPOT ("HOTSPOT/database/").')
    help_info()
    exit(1)

if(os.path.exists(inputs.contigs)==False):
    print(f'Cannot find the input FASTA file. Use "--contigs" to specify the input plasmid contigs')
    help_info()
    exit(1)


#############################################################
##################  Filter short contigs  ###################
#############################################################
rec = []
length_dict = {}
for record in SeqIO.parse(inputs.contigs, 'fasta'):
    length_dict[record.id] = len(record.seq)
    if len(record.seq) > inputs.len:
        rec.append(record)
if(len(rec)==0):
    print(f'No filtered contig left! (with length larger than MINIMUM_LEN)')
    exit(1)
SeqIO.write(rec, f'{out_fn}/filtered_contigs.fa', 'fasta')
pkl.dump(length_dict, open(f'{out_fn}/length_dict', 'wb'))


#############################################################
#####################  Preprocessing  #######################
#############################################################
start = time.time()
threads = inputs.threads
preprocessing(f'{out_fn}/filtered_contigs.fa', db_dir, out_fn, threads)
sentence = assemble_sentences(out_fn, db_dir)
if(len(sentence)==0):
    print(f'The converted sentences are empty! (No available protein)')
    exit(1)


end = time.time()
print(f'The total time for preprocessing is {end-start}s.')
