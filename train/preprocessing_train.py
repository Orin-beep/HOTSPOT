import argparse
import os, csv, sys, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pickle as pkl
import random
from treelib import Tree
import time
start = time.time()


#############################################################
########################  Parameters  #######################
#############################################################
parser = argparse.ArgumentParser(description="""HOTSPOT is a python library for plasmid host prediction.
                                 HOTSPOT is a Transformer-based model and rely on protein-based vocabulary to convert DNA sequences into sentences for prediction.""")
parser.add_argument('--fasta', help='FASTA file of plasmid DNA sequences for training (preferably complete sequences), default: ../training_dataset/plasmids.fasta', default = '../training_dataset/plasmids.fasta')
parser.add_argument('--host_info', help='TSV file containing complete host lineage information from phylum to species for each training plasmid, default: ../training_dataset/host_lineages.tsv', default = '../training_dataset/host_lineages.tsv')
parser.add_argument('--database', help='path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: ../database', type=str, default='../database')
parser.add_argument('--model_path', help='folder to store your customized models, default: models', type=str, default='models')
parser.add_argument('--midfolder', help='folder to store the intermediate files, default: temp', type=str, default='temp')
parser.add_argument('--len', help='minimum length of plasmid DNA sequences, default: 1500', type=int, default=1500)
parser.add_argument('--train_val_list', help='TXT file containing the information of the training/validation sets. The first row should display the list of training plasmids, while the second row should display the list of validation plasmids. Each plasmid in the lists should be separated by a space, default: ../training_dataset/train_val.txt', type=str, default='../training_dataset/train_val.txt')
parser.add_argument('--train_ratio', help='the ratio of the training set size to the total number of input training plasmids. If the train_val_list file is not provided, the training plasmids will be randomly split into train/validation sets using the specified ratio, default: 0.8', type=float, default=0.8)
parser.add_argument('--labels', help="TXT file containing the specified taxonomic labels for different levels. The file should comprise six rows, where each row corresponds to the labels for the phylum, class, order, family, genus, and species levels, respectively. Within each row, the labels should be separated by tabs ('\t'), default: None", type=str, default=None)
parser.add_argument('--num_plasmids', help='minimun number of training plasmids associated with a taxonomic label (used when the labels file is not provided), default: 20', type=int, default=20)
parser.add_argument('--add_frags', help='whether to augment the training set by randomly cutting fragments ranging from 1.5 to 15 kbp (this may slightly slow down the training process, but it will significantly enhance the performance of host prediction), default: True', type=str, default='True')
parser.add_argument('--num_frags', help='maximum number of added fragments from each training plasmid, default: 5', type=int, default=5)
parser.add_argument('--threads', help='number of threads utilized for preprocessing, default: 2', type=str, default='2')
inputs = parser.parse_args()


#############################################################
########################  Help info  ########################
#############################################################
def help_info():
    print('')
    print("""Usage of preprocessing_train.py:
        [--fasta FASTA]    FASTA file of plasmid DNA sequences for training (preferably complete sequences, default: '../training_dataset/plasmids.fasta')
        [--host_info HOST_INFO] TSV file containing complete host lineage information from phylum to species for each training plasmid, default: '../training_dataset/host_lineages.tsv'
        [--database DATABASE]    path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: '../database'
        [--model_path MODEL_PATH]  folder to store your customized models, default: 'models'
        [--midfolder MIDFOLDER]   folder to store the intermediate files, default: 'temp'
        [--len LEN] minimum length of plasmid DNA sequences, default: 1500
        [--train_val_list TRAIN_VAL_LIST]  TXT file containing the information of the training/validation sets. The first row should display the list of training plasmids, while the second row should display the list of validation plasmids. Each plasmid in the lists should be separated by a space, default: '../training_dataset/train_val.txt'
        [--train_ratio TRAIN_RATIO] the ratio of the training set size to the total number of input training plasmids. If the train_val_list file is not provided, the training plasmids will be randomly split into train/validation sets using the specified ratio, default: 0.8
        [--labels LABELS]  TXT file containing the specified taxonomic labels for different levels. The file should comprise six rows, where each row corresponds to the labels for the phylum, class, order, family, genus, and species levels, respectively. Within each row, the labels should be separated by tabs ('\t'), default: None
        [--num_plasmids NUM_PLASMIDS]    minimun number of training plasmids associated with a taxonomic label (used when the labels file is not provided), default: 20
        [--add_frags ADD_FRAGS]   whether to augment the training set by randomly cutting fragments ranging from 1.5 to 15 kbp (this may slightly slow down the training process, but it will significantly enhance the performance of host prediction), default: True
        [--num_frags NUM_FRAGS]   maximum number of added fragments from each training plasmid, default: 5
        [--threads THREADS] number of threads utilized for preprocessing, default: 2
    """)


#############################################################
######################  Check folders  ######################
#############################################################
out_fn = inputs.midfolder
mdl_fn = inputs.model_path
db_fn = inputs.database
for fn in [out_fn, mdl_fn]:
    if not os.path.isdir(fn):
        os.makedirs(fn)


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


#############################################################
###############  Determine taxonomic labels  ################
#############################################################
# read host lineage information
# lineage_dict: lineage dictionary (e.g., {plasmids1: [phylum, class, order, family, genus, species]})
host_info_path = inputs.host_info
lineage_dict = {}
with open(host_info_path, "r") as tsv_file:
    reader = csv.reader(tsv_file, delimiter="\t")
    for row in reader:
        plasmid = row[0]
        if(plasmid not in train_plasmid):
            continue
        lineage_dict[plasmid] = row[1:]
pkl.dump(lineage_dict, open(f'{mdl_fn}/train_lineage.dict', 'wb'))

# check lineage dict and get parent/level info
# parent_dict: e.g., {class1: phylum1} (class1's parent taxon)
# level_dict: e.g., {class1: 1, phylum3: 0, family1: 3}
level_names = {0:'Phylum', 1:'Class', 2:'Order', 3:'Family', 4:'Genus', 5:'Species'}
parent_dict, level_dict = {}, {}
for plasmid in lineage_dict:
    for idx in range(6):
        taxon = lineage_dict[plasmid][idx]
        if(taxon in level_dict and idx!=level_dict[taxon]):
            print(f'Error! Host taxon {taxon} is at two different taxonomic levels ({level_names[level_dict[taxon]]} and {level_names[idx]}) as shown in {host_info_path}.')
            help_info()
            sys.exit()
        level_dict[taxon] = idx
        if(idx==0): # phylum
            parent_dict[taxon] = 'Bacteria'
        else:
            parent = lineage_dict[plasmid][idx-1]
            if(taxon in parent_dict and parent_dict[taxon]!=parent):
                print(f'Error! Host taxon {taxon} has two different parents ({parent} and {parent_dict[taxon]}) as shown in {host_info_path}.')
                help_info()
                sys.exit()
            parent_dict[taxon] = parent
pkl.dump(parent_dict, open(f'{mdl_fn}/parent.dict', 'wb'))

# read the train/validation plasmids, if they are available, or alternatively, randomly split the training plasmids into train/validation sets
# TRAIN: training set
# VAL: validation set
train_val_path = inputs.train_val_list
TRAIN, VAL = set(), set()   # training/validation sets
if(os.path.exists(train_val_path)):
    print('Load the training/validation sets ...')
    f=open(train_val_path)
    ls=f.readlines()
    try:
        TRAIN = set(ls[0].rstrip().split(' '))
        VAL = set(ls[1].rstrip().split(' '))
        f.close()
    except:
        print(f'Error! The train_val_list file should follow the specified format: the first row should list the training plasmids, and the second row should list the validation plasmids. Besides, each plasmid in the lists should be separated by a space.')
        help_info()
        sys.exit()
else:
    train_ratio = inputs.train_ratio
    val_ratio = 1-train_ratio
    print(f'Since the train_val_list file is not provided, the training/validation sets will be randomly split using an {round(train_ratio*100, 2)}% train and {round(val_ratio*100, 2)}% validation ratio.')
    train_size = int(len(train_plasmid)*train_ratio)
    TRAIN = set(random.sample(train_plasmid, train_size))
    VAL = train_plasmid-TRAIN

# write into files
def write_plas_train_val(tmp, out_file):
    f=open(out_file, 'w')
    tmp = sorted(list(tmp))
    for plas in tmp:
        f.write(f'{plas}\n')
    f.close()
write_plas_train_val(TRAIN, f'{mdl_fn}/train.txt')
write_plas_train_val(VAL, f'{mdl_fn}/validation.txt')

# determine taxonomic labels
# labels: label dictionary: {0: {phylum1:0, phylum2:1, ...}}
labels = defaultdict(list)
label_path = inputs.labels
if(label_path!= None and os.path.exists(label_path)):
    print('Load the specified labels ...')
    f=open(label_path)
    ls=f.readlines()
    try:
        for level in range(6):
            labels_tmp = ls[level].rstrip().split('\t')
            labels[level] = sorted(labels_tmp)
            print(f'{level_names[level]}: {len(labels[level])} taxa.')
    except:
        print("Error! The labels file should follow the specified format: the file should comprise six rows, where each row corresponds to the labels for the phylum, class, order, family, genus, and species levels, respectively. Within each row, the labels should be separated by tabs ('\t').")
        help_info()
        sys.exit()
else:
    print(f'Since the labels file is not provided, the labels will be automatically determined based on taxa that have a minimum of {inputs.num_plasmids} associated plasmids in the training set.')
    print('Determine taxonomic labels for training ...')
    label_count = defaultdict(int)
    for plasmid in lineage_dict:
        for taxon in lineage_dict[plasmid]:
            label_count[taxon]+=1
    print(f'Original taxa in {host_info_path}: {len(label_count)}.')

    for taxon in list(label_count):
        if(label_count[taxon]<inputs.num_plasmids):
            del label_count[taxon]
    print(f'After filtering, the taxonomic labels linked to >={inputs.num_plasmids} plasmids: {len(label_count)}.')

    for label in label_count:
        level = level_dict[label]
        labels[level].append(label)
    for level in labels:
        labels[level] = sorted(labels[level])
        print(f'{level_names[level]}: {len(labels[level])} taxa.')

# write in txt and pkl
if not os.path.isdir(f'{mdl_fn}/label_info/'):
    os.makedirs(f'{mdl_fn}/label_info/')
for level in labels:
    tsv_fn = f'{mdl_fn}/label_info/{level_names[level]}.tsv'
    data = [['taxon_name', 'label_id']]

    for idx,taxon in enumerate(labels[level]):
        data.append([taxon, idx]) 

    with open(tsv_fn, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerows(data)

print(f'Detailed label information of each taxonomic level has been written into the folder {mdl_fn}/label_info/.')

# write in a phylogenetic tree with treelib
tree = Tree()
tree.create_node("Bacteria", "Bacteria")
for level in range(6):
    for taxon in labels[level]:
        tree.create_node(taxon, taxon, parent=parent_dict[taxon])
print('The phylogenetic tree of hosts is shown below:')
print(tree)
print(f'The phylogenetic tree of hosts has been saved as {mdl_fn}/host_phylogenetic_tree.pkl.')
pkl.dump(tree, open(f'{mdl_fn}/host_phylogenetic_tree.pkl', 'wb'))

# labels for model training (Labels that lack siblings will be excluded from the set of labels used for Transformer-based models)
labels_mdl = defaultdict(list)
for node in tree.all_nodes():
    taxon = node.identifier
    siblings = tree.siblings(taxon)
    if(len(siblings)!=0):
        level = level_dict[taxon]
        labels_mdl[level].append(taxon)

for level in labels_mdl:    
    labels_mdl[level] = sorted(labels_mdl[level])
    tmp = {}
    for idx,taxon in enumerate(labels_mdl[level]):
        tmp[taxon] = idx
    labels_mdl[level] = tmp
pkl.dump(labels_mdl, open(f'{mdl_fn}/labels.dict', 'wb'))


#############################################################
####################  Data augmentation  ####################
#############################################################
# get plasmids with at least one label
label_set = set()
for level in labels:
    label_set = label_set|set(labels[level])
plasmid_tmp = train_plasmid.copy()
for plasmid in plasmid_tmp:
    if(set(lineage_dict[plasmid])&label_set==set()):
        train_plasmid.remove(plasmid)
rec_dict = {x:y for x,y in rec_dict.items() if x in train_plasmid}
print(f'{len(train_plasmid)} plasmids have at least 1 taxonomic label.')

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# add fragments to rec_dict from train_plasmid
add_frags = str2bool(inputs.add_frags)
if(add_frags==True):
    print('Data augmentation by cutting fragments ...')
    for plasmid in train_plasmid:
        Len = len(rec_dict[plasmid])
        idx=0
        for _ in range(inputs.num_frags):   # default: generate 5 fragments
            rand_len = random.randint(1500, 15000)
            if(rand_len>Len):
                continue
            frag_id = f'{plasmid}+{idx}'
            idx+=1
            init = random.randint(0, Len-rand_len)
            seq = rec_dict[plasmid].seq[init:init+rand_len]
            rec_dict[frag_id] = SeqRecord(Seq(seq), id=frag_id, description='')
    print(f'Data augmentation finished. {len(train_plasmid)} raw plasmids have been augmented to {len(rec_dict)} sequences.')
rec_dict = list(rec_dict.values())
SeqIO.write(rec_dict, f'{out_fn}/train_plasmids.fasta', 'fasta')
del rec_dict


#############################################################
#########################  Prodigal  ########################
#############################################################
print("Running Prodigal ...")
prodigal_cmd = f'prodigal -i {out_fn}/train_plasmids.fasta -a {out_fn}/train_plasmids.faa -f gff -p meta'
_ = subprocess.check_call(prodigal_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


#############################################################
##########################  BLAST  ##########################
#############################################################
thread_num = inputs.threads
print(f"Running alignments with {thread_num} threads ...")
# PC alignment
pc_cmd = f'diamond blastp --threads {thread_num} --sensitive -d {db_fn}/PC_proteins.dmnd -q {out_fn}/train_plasmids.faa -o {out_fn}/resp_pc.tab -k 1'
pc_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resp_pc.tab > {out_fn}/resp_pc.tab.abc"
pc_rm_cmd = f'rm {out_fn}/resp_pc.tab'

# Inc group typing
inc_cmd = f'blastn -query {out_fn}/train_plasmids.fasta -db {db_fn}/database -outfmt 6 -out {out_fn}/resn_inc.tab -num_threads {thread_num} -evalue 1e-10'
inc_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resn_inc.tab > {out_fn}/resn_inc.tab.abc"
inc_rm_cmd = f'rm {out_fn}/resn_inc.tab'

# MOB typing
mob_cmd = f'diamond blastp --threads {thread_num} --sensitive -d {db_fn}/databasemob.dmnd -q {out_fn}/train_plasmids.faa -o {out_fn}/resp_mob.tab -k 1'
mob_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resp_mob.tab > {out_fn}/resp_mob.tab.abc"
mob_rm_cmd = f'rm {out_fn}/resp_mob.tab'
mob_hmm_cmd = f'hmmscan --cpu {thread_num} --incE 0.01 --incdomE 0.01 --domtblout {out_fn}/res_mob.log {db_fn}/MOBfamDB {out_fn}/train_plasmids.faa'

# MPF typing
mpf_cmd = f'diamond blastp --threads {thread_num} --sensitive -d {db_fn}/databasempf.dmnd -q {out_fn}/train_plasmids.faa -o {out_fn}/resp_mpf.tab -k 1'
mpf_awk_cmd = f"awk '{{print $1,$2,$11}}' {out_fn}/resp_mpf.tab > {out_fn}/resp_mpf.tab.abc"
mpf_rm_cmd = f'rm {out_fn}/resp_mpf.tab'

# Run!
cmd_list = [pc_cmd, pc_awk_cmd, pc_rm_cmd, inc_cmd, inc_awk_cmd, inc_rm_cmd, mob_cmd, mob_awk_cmd, mob_rm_cmd, mob_hmm_cmd, mpf_cmd, mpf_awk_cmd, mpf_rm_cmd]
for cmd in cmd_list:
    _ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


#############################################################
####################  Data for training  ####################
#############################################################
# plasmid IDs to be trained at each taxonomic level
level_plas_dict = defaultdict(set)
for plas in train_plasmid:
    for level in labels_mdl:
        if(lineage_dict[plas][level] in labels_mdl[level]):
            level_plas_dict[level].add(plas)
pkl.dump(level_plas_dict, open(f'{mdl_fn}/level_plas.dict', 'wb'))

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
pkl.dump(pc_dict, open(f'{out_fn}/train_feat.dict', 'wb'))
end = time.time()
print(f'Total running time of training preprocessing is {end-start}s.')
