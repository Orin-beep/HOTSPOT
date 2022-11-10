from collections import defaultdict
from bidict import bidict
import csv
import numpy as np
import pickle as pkl
import argparse
import time
import os
import torch
from torch import nn
import torch.utils.data as Data
import sys
sys.path.append('library')
from library.model import hotspot_model


start=time.time()
parser = argparse.ArgumentParser(description="""Introduction: HOTSPOT is a learning-based tool to predict host information from phylum to species for complete plasmids or plasmid contigs assembled from metagenomic data. Its backbone is a phylogenetic tree of the plasmid hosts (bacteria) from phylum to species. By incorporating the state-of-the-art language model, Transformer, in each node’s taxon classifier, the top-down tree search can accurately predict the host taxonomy for the input plasmid contigs. There are totally 115 taxon classifiers, each corresponding to a node with more than one child node. To use HOTSPOT, you only need to input complete plasmids or plasmid contigs assembled from metagenomic data into the program. 
                                                        """)
parser.add_argument('--midfolder', help='folder to store the intermediate files of preprocessing (optional, default temporary_files/)', type=str, default='temporary_files/')
parser.add_argument('--threads', help='number of threads to use (default 8)', type=int, default=8)
parser.add_argument('--mdldir', help='pre-trained model directory (optional)',  default = 'models/')
parser.add_argument('--dbdir', help='database directory (optional, default database/)',  default = 'database/')
parser.add_argument('--out', help='path of the output file (optional, default Result/prediction.tsv)',  type=str, default = 'Result/prediction.tsv')
parser.add_argument('--mode', help='three early stop modes with different uncertainty cutoff estimated with Monte Carlo dropout (MC-dropout). 1: sensitive mode (no early stop used), 2: specific mode (enabling the early stop), 3: accurate mode (enabling the early stop with more stringent uncertainty cutoff, leading to more accurate prediction but returning taxa in higher levels for some inputs). (default: sensitive mode)',  type=int, default = 1)
parser.add_argument('--mcnum', help='the number of the dropout-enabled forward passes to estimate the uncertainty (default: 100, minimum: 10).',  type=int, default = 100)
inputs = parser.parse_args()


#############################################################
######################  Check folders  ######################
#############################################################
def help_info():
    print('')
    print("""The usage of HOTSPOT.py:
            [--midfolder DIR]   Intermediate file folder output by preprocessing.py (default temporary_files/)
            [--mdldir DR]       Path to store the pre-trained models (default models/)
            [--dbdir DR]        Path to store the database (default database/)
            [--out OUT]         Path to store the prediction results (default "Result/prediction.tsv")
            [--threads NUM]     Number of threads to run if 'cpu' is detected ('cuda' not found) (default 8)
            [--mode MOD]        Early stop modes. If 2 or 3 is chosen, the prediction process will slightly slow down.
                                1: sensitive mode (no early stop used)  (default)
                                2: specific mode (enabling the early stop)
                                3: accurate mode (enabling the early stop with more stringent uncertainty cutoff, leading to more accurate prediction but returning taxa in higher levels for some inputs)
                                (default 1)
            [--mcnum MC]        The number of the dropout-enabled forward passes to estimate the uncertainty (default: 100, minimum: 10)
    """)


mdl_dir = inputs.mdldir
if not os.path.exists(mdl_dir):
    print(f'Model directory "{mdl_dir}" missing or unreadable! Please use option "--mdldir" to specify the model path or place the model files under the default path "models/".')
    help_info()
    exit(1)

db_dir = inputs.dbdir
if not os.path.exists(db_dir):
    print(f'Database directory "{db_dir}" missing or unreadable. Please use option "--dbdir" to specify the database path or place the database files under the default path "database/".')
    help_info()
    exit(1)

out_dir = os.path.dirname(inputs.out)
if out_dir != '':
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

out_fn = inputs.midfolder
if(os.path.exists(out_fn)==False):
    print('Cannot find the input intermediate files. Please use "--midfolder" to specify the folder. If you do not have this folder, you can use preprocessing.py to generate it.')
    help_info()
    exit(1)

mcnum = inputs.mcnum
if(mcnum<10):
    mcnum = 10
    print('The specified MC dropout prediction number is smaller than 10. The number will be set to 10.')


#############################################################
######################  Initialization  #####################
#############################################################
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if device == 'cpu':
    print("running with cpu")
    torch.set_num_threads(inputs.threads)

tree = pkl.load(open(f'{db_dir}/tree', 'rb'))
lineage_dict = pkl.load(open(f'{db_dir}/lineage_dict', 'rb'))   # e.g., 'CP057700': ['Escherichia coli', 'Escherichia', 'Enterobacteriaceae', 'Enterobacterales', 'Gammaproteobacteria', 'Proteobacteria', 'Bacteria']
taxa2level = {} # (0, species), (1, genus), (2, family), (3, order), (4, class), (5, phylum)
id2level = {0:'species', 1:'genus', 2:'family', 3:'order', 4:'class', 5:'phylum'}
for t in lineage_dict:
    for i in range(len(lineage_dict[t])):
        taxa2level[lineage_dict[t][i]] = i

mode = inputs.mode
if(mode not in [1,2,3]):
    print(f'The mode input {mode} is not correct. Default 1 (sensitive mode) is used.')
    mode = 1
mode_map = {1: 'sensitive', 2: 'specific', 3: 'accurate'}
blocked = set([])   # early stop contigs
if(mode!=1):
    minimum_var = 0.005
    var_param = pkl.load(open(f'{db_dir}/early_stop_params{mode}', 'rb'))
    for i in var_param:
        var_param[i] = {x:(y if y>=minimum_var else minimum_var) for x,y in var_param[i].items()}


def return_batch(train_sentence, label, flag):
    X_train = torch.from_numpy(train_sentence).to(device)
    y_train = torch.from_numpy(label).to(device)
    train_dataset = Data.TensorDataset(X_train, y_train)
    training_loader = Data.DataLoader(
        dataset=train_dataset,
        batch_size=200,
        shuffle=flag,
        num_workers=0,  # data will be loaded in the main process.
    )
    return training_loader


def enable_dropout(model):
    """ Function to enable the dropout layers during test-time """
    for m in model.modules():
        if m.__class__.__name__.startswith('Dropout'):
            m.train()


#############################################################
########################  Test Data  ########################
#############################################################
X = pkl.load(open(f'{out_fn}/sentence.feat', 'rb'))
id2contig = pkl.load(open(f'{out_fn}/id2contig.dict', 'rb'))    # from 0
id2contig = bidict(id2contig)
length_dict = pkl.load(open(f'{out_fn}/length_dict', 'rb')) # e.g., {'AF318175.1': 1065, 'NZ_CP045289.2': 25142, 'Z22927.1': 3054}
softmax = nn.Softmax(dim=1)
orders = ['phylum', 'class', 'order', 'family', 'genus', 'species'] # the output format
models = os.listdir(mdl_dir)

res = defaultdict(list) # prediction results
preds_tmp = defaultdict(list)   # contigs corresponding to currect model
preds_tmp['Bacteria'] = list(range(len(id2contig))) # all plasmids initially
pending = ['Bacteria']
print(f'Run HOTSPOT in {mode_map[mode]} mode ...')


#############################################################
########################  Prediction  #######################
#############################################################
while(pending!=[]):
    model_name = pending[0]
    #print('_______________________________________')
    print('\n######################################################################')
    print(f'Identifying classes in {model_name} ...')

    # if only one child, no model needed
    if(f'{model_name}.pth' not in models):
        del pending[0]
        if(len(tree.children(model_name))==1):
            for i in preds_tmp[model_name]:
                if(id2contig[i] not in blocked):
                    res[id2contig[i]].append(tree.children(model_name)[0].identifier)
            taxa = tree.children(model_name)[0].identifier
            pending.append(taxa)
            preds_tmp[taxa] = preds_tmp[model_name]
            print(id2level[taxa2level[taxa]])
            print(f'{len(preds_tmp[model_name])} plasmid contigs classified into {id2level[taxa2level[taxa]]} {taxa}')
        else:
            print(f'Achieve leaf node "{id2level[taxa2level[model_name]]} {model_name}"!')
        continue

    # run HOTSPOT model
    contigs = []
    idxx = []
    test = set(preds_tmp[model_name])
    for i in range(len(id2contig)):
        if(i in test):
            contig = id2contig[i]
            idxx.append(i)
            contigs.append(contig)
    X_test = X[idxx]
    y_test = np.array([0]*X_test.shape[0])
    print(f'{len(X_test)} plasmid contigs to be predicted ...')
    mapping = pkl.load(open(f'{db_dir}/Labels_mapping/{model_name}', 'rb'))

    # load model
    model = hotspot_model(device=device,dropout=0.5,output_num = len(mapping)).to(device)
    model.load_state_dict(torch.load(f'{mdl_dir}/{model_name}.pth',map_location=device))
    _ = model.eval()
    test_loader = return_batch(X_test, y_test, flag = False)
    with torch.no_grad():
        preds = []
        for step, (batch_x, batch_y) in enumerate(test_loader):
            sentence = batch_x.int()
            logit, _ = model(sentence)
            logit = softmax(logit)
            pred = [torch.argmax(item) for item in logit]
            preds+=pred
        idx=0

    # dropout enabled
    if(mode != 1):
        print(f'Computing Monte Carlo Dropout method {mcnum} times ...')
        dropout_predictions = np.empty((0, X_test.shape[0], len(mapping)))
        for i in range(mcnum):    # {mcnum} MC dropout
            predictions =np.empty((0, len(mapping)))
            model.eval()
            enable_dropout(model)
            for step, (batch_x, batch_y) in enumerate(test_loader):
                sentence = batch_x.int()
                with torch.no_grad():
                    logit, _ = model(sentence)
                    logit = softmax(logit)
                predictions = np.vstack((predictions, logit.cpu().numpy()))
            dropout_predictions = np.vstack((dropout_predictions,predictions[np.newaxis, :, :]))
        variance = np.var(dropout_predictions, axis=0) # shape (n_samples, n_classes), e.g. (4, 10) 4 plasmid contigs and 10 classes

    # generate results
    count_num = defaultdict(int)
    block_num = 0
    idx=0
    for i in preds:
        r = mapping.inv[i.item()]
        contig = contigs[idx]

        # judge variance
        if(mode!=1):
            L = length_dict[contig]
            if(L<2000):
                num=1500
            elif(L<4000):
                num=3000
            elif(L<7500):
                num=5000
            else:
                num=10000
            var_cutoff = var_param[model_name][num]
            if(variance[idx][i.item()]>=var_cutoff):
                blocked.add(contig)

        if(contig not in blocked):
            res[contig].append(r)
            count_num[r]+=1
            preds_tmp[r].append(idxx[idx])
            if(r not in pending):
                pending.append(r)
        else:
            block_num+=1
        idx+=1

    for i in count_num:
        print(f'{count_num[i]} plasmid contigs classified into {id2level[taxa2level[i]]} {i}.')
    if(block_num!=0):
        print(f'{block_num} plasmid contigs stopping prediction with MC Dropout!)')

    del pending[0]


#############################################################
#########################  Results  #########################
#############################################################
tsv_w = csv.writer(open(inputs.out, 'w'), delimiter='\t')
tsv_w.writerow(['Contig', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
for i in res:
    tmp = ['-']*7
    tmp[0] = i
    for j in range(len(res[i])):
        tmp[j+1] = res[i][j]
    tsv_w.writerow(tmp)

end=time.time()
print(f"""Plasmid host prediction finished!
The prediction results are saved in {inputs.out} (total running time: {end-start}s).""")
