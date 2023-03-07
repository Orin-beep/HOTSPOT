import pickle as pkl
import numpy as np
import pandas as pd
from collections import defaultdict
import os


def read_replicon(blastn_file):   # read incompatibility group features
    inc_list = ['IncA/C', 'IncY', 'IncP', 'IncHI1', 'Inc4', 'IncHI2', 'IncI/B/O/K/Z', 'IncT', 'IncQ', 'Inc11', 'Inc13', 'IncN', 'IncFIC', 'FII', 'IncU', 'IncW', 'IncX', 'IncL/M', 'IncR', 'FIA', 'IncFIB', 'Inc18']   # 22 inc groups
    incIdx_dict = {x:y for y,x in enumerate(inc_list)}
    f=open(blastn_file)
    ls=f.readlines()
    contig2rep = defaultdict(dict)
    for l in ls:
        d = l.rstrip().split(' ')
        rep = d[1].split('_')[0]
        contig2rep[d[0]][incIdx_dict[rep]] = 1
    return contig2rep


def read_mpf_mob(blastp_file, hmmer_file):  # read mating pair formation (MPF), and relaxase type (MOB)
    # do not consider unknown mob/mpf type
    blastp_idx = {'MOBB':4, 'MOBQ':7, 'MOBP':8, 'MOBM':5, 'MOBF':1, 'MOBT':2, 'MOBC':3, 'MOBH':6, 'MOBV':9, 'MPF_G':10, 'MPF_T':11, 'MPF_F':12, 'MPF_I':13}
    hmm_idx = {'profile_MOBF':1, 'profile_MOBT':2, 'T4SS_MOBC':3, 'T4SS_MOBB':4, 'profile_MOBM':5, 'T4SS_MOBH':6, 'T4SS_MOBQ':7, 'T4SS_MOBP3':8, 'T4SS_MOBP1':8, 'T4SS_MOBP2':8, 'T4SS_MOBV':9}
    res_dict = defaultdict(dict)

    # read MOBscan results
    f=open(hmmer_file)
    ls=f.readlines()
    for l in ls:
        if(l[0]=='#'):
            continue
        d=l.rstrip().split(' ')
        d=[x for x in d if x!='']   # 23 columns totally
        mob=hmm_idx[d[0]]
        query=d[3]
        Evalue = 0  # hmmscan results are of the highest priority
        res_dict[query] = (mob, Evalue) if query not in res_dict or Evalue<res_dict[query][1] else res_dict[query]

    # read blastp results
    f=open(blastp_file)
    ls=f.readlines()
    for l in ls:
        d=l.rstrip().split(' ')
        query=d[0]
        res=d[1].split('|')[1]
        if(res not in blastp_idx):
            continue
        res = blastp_idx[res]
        Evalue = float(d[2])
        res_dict[query] = (res, Evalue) if query not in res_dict or Evalue<res_dict[query][1] else res_dict[query]

    contig2mob = defaultdict(list)
    for i in res_dict:
        contig = i.rsplit('_', 1)[0]
        idx = int(i.rsplit('_', 1)[1])
        contig2mob[contig].append((idx, res_dict[i][0]))
    for contig in contig2mob:   # sorting based on protein positions
        contig2mob[contig] = sorted(contig2mob[contig], key=lambda tup: tup[0])

    return contig2mob


def assemble_sentence(dbPath, mclBlastp_file, mclLength, out_fn, contig2mob, mob_mpfLength, contig2rep):
    proteins_df = pd.read_csv(dbPath+'/proteins.csv')
    proteins_df.dropna(axis=0, how='any', inplace=True)
    pc2wordsid = {pc: idx for idx, pc in enumerate(sorted(set(proteins_df['cluster'].values)))} # MCL protein cluster index dictionary
    max_num = max(list(pc2wordsid.values()))
    protein2pc = {protein: pc for protein, pc in zip(proteins_df['protein_id'].values, proteins_df['cluster'].values)}
    blast_df = pd.read_csv(mclBlastp_file, sep=' ', names=['query', 'ref', 'evalue'])

    contig2pcs = defaultdict(list)
    for query, ref, evalue in zip(blast_df['query'].values, blast_df['ref'].values, blast_df['evalue'].values):
        contig = query.rsplit('_', 1)[0]
        idx = int(query.rsplit('_', 1)[1])
        try:
            pc = pc2wordsid[protein2pc[ref]]
        except:
            continue
        contig2pcs[contig].append((idx, pc, evalue))

    # sorting based on protein positions
    for contig in contig2pcs:
        contig2pcs[contig] = sorted(contig2pcs[contig], key=lambda tup: tup[0])

    # Contigs2sentence
    contig2id = {contig:idx for idx, contig in enumerate(contig2pcs.keys())}
    id2contig = {idx:contig for idx, contig in enumerate(contig2pcs.keys())}
    sentence = np.zeros((len(contig2id.keys()), mclLength))
    for row in range(sentence.shape[0]):
        contig = id2contig[row]
        pcs = contig2pcs[contig]
        for col in range(len(pcs)):
            try:
                _, sentence[row][col], _ = pcs[col]
                sentence[row][col] += 1 # minimal number is 1
            except:
                break
    sentence2 = np.zeros((len(id2contig), mob_mpfLength))
    for row in range(sentence2.shape[0]):
        contig = id2contig[row]
        mob = contig2mob[contig]
        for col in range(len(mob)):
            try:
                _, sentence2[row][col] = mob[col]
            except:
                break
    sentence3 = np.zeros((len(id2contig), 22))
    for row in range(sentence3.shape[0]):
        contig = id2contig[row]
        rep = contig2rep[contig]
        for col in rep:
            try:
                sentence3[row][col] = rep[col]
            except:
                break

    sentence = np.hstack((sentence, sentence2, sentence3))
    print('The encoded sentence matrix:')
    print(sentence)
    print("The matrix's shape is:", end = ' ')
    print(sentence.shape)
    num_pcs = len(set(pc2wordsid.values()))

    pkl.dump(sentence, open(f'{out_fn}/sentence.feat', 'wb'))
    pkl.dump(id2contig, open(f'{out_fn}/id2contig.dict', 'wb'))
    return sentence


def assemble_sentences(fn, mcldb_path):
    contig2rep = read_replicon(f'{fn}/resN.tab.abc')
    contig2mob = read_mpf_mob(f'{fn}/resM.tab.abc', f'{fn}/res.log')
    sentence = assemble_sentence(mcldb_path, f'{fn}/resp.tab.abc', 400, fn, contig2mob, 50, contig2rep)
    return sentence


def preprocessing(fa, database, out_fn, t):
    # blastn against inc group databases
    os.system(f'blastn -query {fa} -db {database}/database -outfmt 6 -out {out_fn}/resN.tab -num_threads {t}')
    os.system(f"awk '{{print $1,$2,$11}}' {out_fn}/resN.tab > {out_fn}/resN.tab.abc")
    os.system(f'rm {out_fn}/resN.tab')

    # blastp and hmmscan against mob/mpf databases
    os.system(f'prodigal -i {fa} -a {out_fn}/plasmids.faa -f gff -p meta')
    os.system(f'diamond blastp --threads {t} --sensitive -d {database}/databasemob.dmnd -q {out_fn}/plasmids.faa -o {out_fn}/resmob.tab -k 1')
    os.system(f"awk '{{print $1,$2,$11}}' {out_fn}/resmob.tab > {out_fn}/resmob.tab.abc")
    os.system(f'rm {out_fn}/resmob.tab')

    os.system(f'diamond blastp --threads {t} --sensitive -d {database}/databasempf.dmnd -q {out_fn}/plasmids.faa -o {out_fn}/resmpf.tab -k 1')
    os.system(f"awk '{{print $1,$2,$11}}' {out_fn}/resmpf.tab > {out_fn}/resmpf.tab.abc")
    os.system(f'rm {out_fn}/resmpf.tab')

    os.system(f'cat {out_fn}/resmpf.tab.abc {out_fn}/resmob.tab.abc > {out_fn}/resM.tab.abc')
    os.system(f'rm {out_fn}/resmpf.tab.abc {out_fn}/resmob.tab.abc')
    os.system(f'hmmscan --cpu {t} --incE 0.01 --incdomE 0.01 --domtblout {out_fn}/res.log {database}/MOBfamDB {out_fn}/plasmids.faa')
    os.system(f'diamond blastp --threads {t} --sensitive -d {database}/database.dmnd -q {out_fn}/plasmids.faa -o {out_fn}/resp.tab -k 1')
    os.system(f"awk '{{print $1,$2,$11}}' {out_fn}/resp.tab > {out_fn}/resp.tab.abc")
    os.system(f'rm {out_fn}/resp.tab {out_fn}/plasmids.faa')
