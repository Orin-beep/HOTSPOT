import argparse
import os, sys, subprocess
from Bio import SeqIO
from collections import defaultdict
import pickle as pkl
import time
start = time.time()


#############################################################
########################  Parameters  #######################
#############################################################
parser = argparse.ArgumentParser(description="""HOTSPOT is a python library for plasmid host prediction.
                                 HOTSPOT is a Transformer-based model and rely on protein-based vocabulary to convert DNA sequences into sentences for prediction.""")
parser.add_argument('--fasta', help='FASTA file of the plasmid DNA sequences to be predicted (either complete sequences or contigs), default: example_plasmids/multiple_plasmids.fasta', default = 'example_plasmids/multiple_plasmids.fasta')
parser.add_argument('--database', help='path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: database', type=str, default='database')
parser.add_argument('--model_path', help='path of the folder storing the downloaded or your customized models, default: models', type=str, default='models')
parser.add_argument('--midfolder', help='folder to store the intermediate files for prediction, default: temp', type=str, default='temp')
parser.add_argument('--len', help='minimum length of plasmid DNA sequences, default: 1500', type=int, default=1500)
parser.add_argument('--threads', help='number of threads utilized for preprocessing, default: 2', type=str, default='2')
inputs = parser.parse_args()


#############################################################
########################  Help info  ########################
#############################################################
def help_info():
    print('')
    print("""Usage of preprocessing.py:
        [--fasta FASTA] FASTA file of the plasmid DNA sequences to be predicted (either complete sequences or contigs), default: multiple_plasmids.fasta
        [--database DATABASE]   path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: database
        [--model_path MODEL_PATH]   path of the folder storing the downloaded or your customized models, default: models
        [--midfolder]   folder to store the intermediate files for prediction, default: temp
        [--len LEN] minimum length of plasmid DNA sequences, default: 1500
        [--threads THREADS] number of threads utilized for preprocessing, default: 2
    """)


#############################################################
######################  Check folders  ######################
#############################################################
out_fn = inputs.midfolder
if not os.path.isdir(out_fn):
    os.makedirs(out_fn)
mdl_fn = inputs.model_path
db_fn = inputs.database

if not os.path.exists(db_fn):
    print(f"Error! The database folder '{db_fn}' is unavailable. Please use the option '--database' to indicate the directory of the downloaded database folder.")
    help_info()
    sys.exit()

if not os.path.exists(mdl_fn):
    print(f"Error! The model folder '{mdl_fn}' is unavailable. Please use the option '--model_path' to indicate the folder of the downloaded or your customized models.")
    help_info()
    sys.exit()


#############################################################
#################  Filter short sequences  ##################
#############################################################
# rec_dict: Seq-record dictionary; 
# train_plasmid: set of all the input plasmid IDs.
rec_dict, train_plasmid = {}, set()
for s in SeqIO.parse(inputs.fasta, 'fasta'):
    if(len(s.seq)>=inputs.len):
        rec_dict[s.id] = s
        train_plasmid.add(s.id)
print(f'{len(train_plasmid)} plasmids have length >={inputs.len} bp.')
if(len(train_plasmid)==0):
    print(f'Error! No input plasmid DNA sequences exceed a length of {inputs.len} bp.')
    help_info()
    sys.exit()

# write into FASTA file in temp folder
rec_dict = list(rec_dict.values())
SeqIO.write(rec_dict, f'{out_fn}/test_plasmids.fasta', 'fasta')


#############################################################
#########################  Prodigal  ########################
#############################################################
print("Running Prodigal ...")
prodigal_cmd = f'prodigal -i {out_fn}/test_plasmids.fasta -a {out_fn}/test_plasmids.faa -f gff -p meta'
_ = subprocess.check_call(prodigal_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


#############################################################
##########################  BLAST  ##########################
#############################################################
thread_num = inputs.threads
print(f"Running alignments with {thread_num} threads ...")
# PC alignment
pc_cmd = f'diamond blastp --threads {thread_num} --sensitive -d {db_fn}/PC_proteins.dmnd -q {out_fn}/test_plasmids.faa -o {out_fn}/resp_pc.tab -k 1'
pc_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resp_pc.tab > {out_fn}/resp_pc.tab.abc"
pc_rm_cmd = f'rm {out_fn}/resp_pc.tab'

# Inc group typing
inc_cmd = f'blastn -query {out_fn}/test_plasmids.fasta -db {db_fn}/database -outfmt 6 -out {out_fn}/resn_inc.tab -num_threads {thread_num} -evalue 1e-10'
inc_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resn_inc.tab > {out_fn}/resn_inc.tab.abc"
inc_rm_cmd = f'rm {out_fn}/resn_inc.tab'

# MOB typing
mob_cmd = f'diamond blastp --threads {thread_num} --sensitive -d {db_fn}/databasemob.dmnd -q {out_fn}/test_plasmids.faa -o {out_fn}/resp_mob.tab -k 1'
mob_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resp_mob.tab > {out_fn}/resp_mob.tab.abc"
mob_rm_cmd = f'rm {out_fn}/resp_mob.tab'
mob_hmm_cmd = f'hmmscan --cpu {thread_num} --incE 0.01 --incdomE 0.01 --domtblout {out_fn}/res_mob.log {db_fn}/MOBfamDB {out_fn}/test_plasmids.faa'

# MPF typing
mpf_cmd = f'diamond blastp --threads {thread_num} --sensitive -d {db_fn}/databasempf.dmnd -q {out_fn}/test_plasmids.faa -o {out_fn}/resp_mpf.tab -k 1'
mpf_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resp_mpf.tab > {out_fn}/resp_mpf.tab.abc"
mpf_rm_cmd = f'rm {out_fn}/resp_mpf.tab'

# Run!
cmd_list = [pc_cmd, pc_awk_cmd, pc_rm_cmd, inc_cmd, inc_awk_cmd, inc_rm_cmd, mob_cmd, mob_awk_cmd, mob_rm_cmd, mob_hmm_cmd, mpf_cmd, mpf_awk_cmd, mpf_rm_cmd]
for cmd in cmd_list:
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


#############################################################
###################  Data for prediction  ###################
#############################################################
# read PC
# pc_dict, plas: [55, 604, 9046, 7, ...] with length of 400
f=open(f'{out_fn}/resp_pc.tab.abc')
ls=f.readlines()
pc_dict = defaultdict(dict)
for l in ls:
    d = l.split(' ')
    query = d[0]
    target = d[1]
    plas = query[:query.rfind('_')]
    idx = int(query[query.rfind('_')+1:])
    pc_id = int(target[target.find('_')+1:target.rfind('_')])+1 # from 1, 0 is padding token
    pc_dict[plas][idx] = pc_id
for plas in pc_dict:
    tmp = sorted(pc_dict[plas].items())
    tmp = [i[1] for i in tmp]
    tmp = tmp[:400]
    tmp = tmp + [0]*(400-len(tmp))
    pc_dict[plas] = tmp

# if there is no any PC aligned
if(len(pc_dict)==0):
    print(f'Error! No PC token aligned for the input plasmid DNA sequences.')
    help_info()
    sys.exit()

# read Inc groups
inc_list = ['IncA/C', 'IncY', 'IncP', 'IncHI1', 'Inc4', 'IncHI2', 'IncI/B/O/K/Z', 'IncT', 'IncQ', 'Inc11', 'Inc13', 'IncN', 'IncFIC', 'FII', 'IncU', 'IncW', 'IncX', 'IncL/M', 'IncR', 'FIA', 'IncFIB', 'Inc18']
INC_dict = {y:x for x,y in enumerate(inc_list)}
f = open(f'{out_fn}/resn_inc.tab.abc')
ls=f.readlines()
inc_dict = {}
for l in ls:
    d = l.split(' ')
    plas = d[0]
    target = d[1]
    inc = target[:target.rfind('_')]
    inc_idx = INC_dict[inc]
    try:
        inc_dict[plas][inc_idx] = 1
    except:
        inc_dict[plas] = [0]*len(inc_list)
        inc_dict[plas][inc_idx] = 1

# read MOB/MPF from Diamond
mob_mpf_dict = {'MOBB':4, 'MOBQ':7, 'MOBP':8, 'MOBM':5, 'MOBF':1, 'MOBT':2, 'MOBC':3, 'MOBH':6, 'MOBV':9, 'MPF_G':10, 'MPF_T':11, 'MPF_F':12, 'MPF_I':13}
f=open(f'{out_fn}/resp_mob.tab.abc')
ls=f.readlines()
f.close()
f=open(f'{out_fn}/resp_mpf.tab.abc')
ls = ls + f.readlines()
f.close()
M_dict = defaultdict(dict)
for l in ls:
    d = l.rstrip().split(' ')
    query = d[0]
    plas = query[:query.rfind('_')]
    idx = int(query[query.rfind('_')+1:])
    target = d[1]
    evalue = float(d[2])
    M = target.split('|')[-1]
    if(M not in mob_mpf_dict):
        continue
    M_idx = mob_mpf_dict[M]
    M_dict[plas][idx] = [M_idx, evalue]

# read MOB from MOBscan
hmm_dict = {'profile_MOBF':1, 'profile_MOBT':2, 'T4SS_MOBC':3, 'T4SS_MOBB':4, 'profile_MOBM':5, 'T4SS_MOBH':6, 'T4SS_MOBQ':7, 'T4SS_MOBP3':8, 'T4SS_MOBP1':8, 'T4SS_MOBP2':8, 'T4SS_MOBV':9}
f=open(f'{out_fn}/res_mob.log')
ls=f.readlines()
for l in ls:
    if(l[0]=='#'):
        continue
    d=l.rstrip().split(' ')
    d=[x for x in d if x!='']
    mob=hmm_dict[d[0]]
    query=d[3]
    plas = query[:query.rfind('_')]
    idx = int(query[query.rfind('_')+1:])
    evalue = float(d[6])
    if(idx not in M_dict[plas]):
        M_dict[plas][idx] = [mob, evalue]
        continue
    if(evalue<M_dict[plas][idx][1]):
        M_dict[plas][idx] = [mob, evalue]
for plas in M_dict:
    for idx in M_dict[plas]:
        M_dict[plas][idx] = M_dict[plas][idx][0]

for plas in M_dict:
    tmp = sorted(M_dict[plas].items())
    tmp = [i[1] for i in tmp]
    tmp = tmp[:50]
    tmp = tmp + [0]*(50-len(tmp))
    M_dict[plas] = tmp

# combine PC, MOB/MPF, and Inc
for plas in pc_dict:
    if(plas not in M_dict):
        pc_dict[plas] += [0]*50
    else:
        pc_dict[plas] += M_dict[plas]
    if(plas not in inc_dict):
        pc_dict[plas] += [0]*len(inc_list)
    else:
        pc_dict[plas] += inc_dict[plas]
    pc_dict[plas] = pc_dict[plas]
pkl.dump(pc_dict, open(f'{out_fn}/test_feat.dict', 'wb'))
end = time.time()
print(f'Total running time of the prediction preprocessing is {end-start}s.')
