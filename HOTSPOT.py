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
parser = argparse.ArgumentParser(description="HOTSPOT is a learning-based tool designed for plasmid host prediction. Its backbone is a phylogenetic tree of plasmids' hosts (bacteria) from phylum to species. The top-down tree search can accurately predict the hosts' taxonomic labels by incorporating the state-of-the-art language model, Transformer, in each node’s taxon classifier. To use HOTSPOT, you only need to input plasmid sequences (complete or segmented) into the program.")
parser.add_argument('--midfolder', help='folder to store the intermediate files from preprocessing (used as the inputs of HOTSPOT.py, default temporary_files/)', type=str, default='temporary_files/')
parser.add_argument('--threads', help="number of threads to use if 'cpu' is detected ('cuda' not found, default 8)", type=int, default=8)
parser.add_argument('--mdldir', help="pre-trained models' directory (default models/)",  default = 'models/')
parser.add_argument('--dbdir', help='database directory (default database/)',  default = 'database/')
parser.add_argument('--out', help='path to store the output files (default Results/)',  type=str, default = 'Results/')
parser.add_argument('--accurate', help='''If this parameter is 1, the MC-dropout based early stop mechanism will be activated with two sets of uncertainty cutoffs, and the prediction will cost more time. 1)sensitive mode (the default mode without early stop, output: 'Results/host_lineage.tsv'). 2)specific mode (enabling the early stop, output: 'Results/host_lineage_specific.tsv'). 3)accurate mode (enabling the early stop with more stringent uncertainty cutoff, leading to more accurate prediction but returning taxa in higher levels for some inputs, output: 'Results/host_lineage_accurate.tsv'). default: 0.'''
        , type=int, default = 0)
parser.add_argument('--mcnum', help='''the number of the dropout-enabled forward passes to estimate the prediction uncertainty (works when '--accurate 1' is chosen, default: 100, minimum: 10).''',  type=int, default = 100)
inputs = parser.parse_args()


#############################################################
######################  Check folders  ######################
#############################################################
def help_info():
    print('')
    print("""The usage of HOTSPOT.py:
            [--midfolder DIR]   Folder to store the intermediate files from preprocessing (used as the inputs of HOTSPOT.py, default temporary_files/)
            [--mdldir DR]       Pre-trained models' directory (default models/)
            [--dbdir DR]        Database directory (default database/)
            [--out OUT]         Path to store the output files (default Results/)
            [--threads NUM]     Number of threads to use if 'cpu' is detected ('cuda' not found, default 8)
            [--accurate ACC]    If this parameter is 1, the MC-dropout based early stop mechanism will be activated with two sets of uncertainty cutoffs, and the prediction will cost more time.
                                1. sensitive mode (the default mode without early stop, output: 'Results/host_lineage.tsv')
                                2. specific mode (enabling the early stop, output: 'Results/host_lineage_specific.tsv')
                                3. accurate mode (enabling the early stop with more stringent uncertainty cutoff, leading to more accurate prediction but returning taxa in higher levels for some inputs, output: 'Results/host_lineage_accurate.tsv')
                                (default 0)
            [--mcnum MC]        The number of the dropout-enabled forward passes to estimate the prediction uncertainty (works when '--accurate 1' is chosen, default: 100, minimum: 10)
    """)


mdl_dir = inputs.mdldir
if not os.path.exists(mdl_dir):
    print(f'''Model directory "{mdl_dir}" unavailable! Please use the option "--mdldir" to specify the model path or place the model folder under HOTSPOT's main directory ("HOTSPOT/models/").''')
    help_info()
    exit(1)

db_dir = inputs.dbdir
if not os.path.exists(db_dir):
    print(f'''Database directory "{db_dir}" unavailable. Please use the option "--dbdir" to specify the database path or place the database folder under HOTSPOT's  main directory ("HOTSPOT/database/").''')
    help_info()
    exit(1)

out_dir = inputs.out
if(os.path.exists(out_dir)==False):
    os.system(f'mkdir {out_dir}')

out_fn = inputs.midfolder
if(os.path.exists(out_fn)==False):
    print('Cannot find the input intermediate files. Please use "--midfolder" to specify the folder. If you do not have this folder, you can use preprocessing.py to generate it.')
    help_info()
    exit(1)

early_stop = inputs.accurate
if(early_stop==1):
    mcnum = inputs.mcnum
    if(mcnum<10):
        mcnum = 10
        print('The specified number of the dropout-enabled forward passes is smaller than 10. The number is set to 10.')
    blocked2, blocked3 = set(), set()
    res2 = defaultdict(list)    # specific
    res3 = defaultdict(list)    # accurate
    print(f'Run HOTSPOT with early stop ...')
else:
    print(f'Run HOTSPOT without early stop ...')


#############################################################
######################  Initialization  #####################
#############################################################
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if device == torch.device('cpu'):
    print("Running with cpu ...")
    torch.set_num_threads(inputs.threads)

tree = pkl.load(open(f'{db_dir}/tree.pkl', 'rb'))
lineage_dict = pkl.load(open(f'{db_dir}/lineage.dict', 'rb'))


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
softmax = nn.Softmax(dim=1)
models = os.listdir(mdl_dir)

res1 = defaultdict(list) # prediction results
preds_tmp = defaultdict(list)   # contigs corresponding to currect model
preds_tmp['Bacteria'] = list(range(len(id2contig))) # all plasmids initially
pending = ['Bacteria']


#############################################################
########################  Prediction  #######################
#############################################################
while(pending!=[]):
    model_name = pending[0]
    print('\n######################################################################')
    print(f'Identifying classes in {model_name} ...')
    children = tree.children(model_name)
    children = sorted([x.identifier for x in children])

    # if only one child, no model needed
    if(len(children)<2): 
        del pending[0]
        if(len(children)==1):
            for i in preds_tmp[model_name]:
                res1[id2contig[i]].append(children[0]+'*')
                if(early_stop==1):
                    if(id2contig[i] not in blocked2):
                        res2[id2contig[i]].append(children[0]+'*')
                    if(id2contig[i] not in blocked3):
                        res3[id2contig[i]].append(children[0]+'*')
            child = children[0]
            pending.append(child)
            preds_tmp[child] = preds_tmp[model_name]
        continue

    # run HOTSPOT model
    label_dict = pkl.load(open(f'{db_dir}/label_dict/{model_name}.pkl', 'rb'))
    label_dict = bidict(label_dict)

    X_idx = sorted(preds_tmp[model_name])
    contigs = []
    for i in X_idx:
        contig = id2contig[i]
        contigs.append(contig)
    X_test = X[X_idx]
    y_test = np.array([0]*X_test.shape[0])
    print(f'{len(X_test)} plasmid contigs to be predicted ...')

    # load model
    model = hotspot_model(device=device,dropout=0.5,output_num = len(label_dict)).to(device)
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

    # dropout enabled
    if(early_stop==1):
        print(f'Computing Monte Carlo Dropout for {mcnum} times ...')
        dropout_predictions = np.empty((0, X_test.shape[0], len(label_dict)))
        for i in range(mcnum):    # {mcnum} MC dropout
            predictions =np.empty((0, len(label_dict)))
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
        proba = np.mean(dropout_predictions, axis=0, keepdims=False)

    # summarize results
    idx=0
    for i in preds:
        r = label_dict.inv[i.item()]
        contig = contigs[idx]
        res1[contig].append(r)
        preds_tmp[r].append(X_idx[idx])
        if(r not in pending):
            pending.append(r)
        if(early_stop==1):
            prob_mean = proba[idx][i.item()]
            VAR_ratio = variance[idx][i.item()]/prob_mean
            if(VAR_ratio>=0.4 or prob_mean<0.5):
                blocked2.add(contig)
            if(VAR_ratio>=0.3 or prob_mean<0.75):
                blocked3.add(contig)
            if(contig not in blocked2):
                res2[contig].append(r)
            if(contig not in blocked3):
                res3[contig].append(r)
        idx+=1
    del pending[0]


#############################################################
#########################  Results  #########################
#############################################################
tsv_w = csv.writer(open(f'{out_dir}/host_lineage.tsv', 'w'), delimiter='\t')
tsv_w.writerow(['Contig', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
for i in res1:
    tmp = ['-']*7
    tmp[0] = i
    for j in range(len(res1[i])):
        tmp[j+1] = res1[i][j]
    tsv_w.writerow(tmp)
if(early_stop==1):
    tsv_w2 = csv.writer(open(f'{out_dir}/host_lineage_specific.tsv', 'w'), delimiter='\t')
    tsv_w3 = csv.writer(open(f'{out_dir}/host_lineage_accurate.tsv', 'w'), delimiter='\t')
    tsv_w2.writerow(['Contig', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    tsv_w3.writerow(['Contig', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    for i in res2:
        tmp = ['-']*7
        tmp[0] = i
        for j in range(len(res2[i])):
            tmp[j+1] = res2[i][j]
        tsv_w2.writerow(tmp)
    for i in res3:
        tmp = ['-']*7
        tmp[0] = i
        for j in range(len(res3[i])):
            tmp[j+1] = res3[i][j]
        tsv_w3.writerow(tmp)


end=time.time()
print(f"""Plasmid host prediction completed!
The prediction results are saved in {inputs.out}.
The total running time for HOTSPOT.py is {end-start}s.""")
