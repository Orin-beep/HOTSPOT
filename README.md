# HOTSPOT: Hierarchical hOst predicTion for aSsembled Plasmid cOntigs with Transformer

HOTSPOT's backbone is a phylogenetic tree of plasmids' hosts (bacteria) from phylum to species. The top-down tree search can accurately predict the hosts' taxonomic labels by incorporating the state-of-the-art language model, Transformer, in each node’s taxon classifier. __To use HOTSPOT, you only need to input plasmid sequences (complete or segmented) into the program.__


# Required Dependencies
* Python 3.x
* Numpy
* bidict
* Pandas
* Pytorch>1.8.0
* [Diamond](https://github.com/bbuchfink/diamond)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Biopython](https://biopython.org/)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [HMMER](http://hmmer.org/)
* [treelib](https://pypi.org/project/treelib/)

If you want to use the gpu to accelerate the program:
* cuda
* Pytorch-gpu
* For cpu version pytorch: ```conda install pytorch torchvision torchaudio cpuonly -c pytorch```
* For gpu version pytorch: Search [pytorch](https://pytorch.org/) to find the correct cuda version according to your computer


## Prepare the environment
After cloning this respository (```git clone https://github.com/Orin-beep/HOTSPOT```), you can use anaconda to install the environment.yaml. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: 
```
conda env create -f environment.yaml -n hotspot
conda activate hotspot
```


## Prepare the database (from Google Drive)
To prepare the database, you can simply use this bash script: 
```
sh prepare_db.sh
```
To be sure that the downloaded `database/` folder is placed in the HOTSPOT's main directory. Otherwise, you have to use option `--dbdir` to specify the database path when running `preprocessing.py` and `HOTSPOT.py`.


## Prepare the pre-trained models (from Google Drive)
If you don't need the host identification function at the species level, you may only download the pre-trained models from phylum to genus (25.6G):
```
sh prepare_light_mdl.sh
```

Or if you have sufficient disk storage, you can download the complete models (43.9G):
```
sh prepare_mdl.sh
```

Similar to the database preparation, if the model path is not under HOTSPOT's main directory, you need to use option `--mdldir` to specify the model path when running `HOTSPOT.py`.

The pre-trained models are a little big because the phylogenetic includes 10, 19, 46, 94, 180, and 272 nodes for the six taxonomic ranks from phylum to species, respectively. We will try to optimize the model sizes in the next version of HOTSPOT.


## Download the pre-trained models manually
If the `prepare_light_mdl.sh` or `prepare_mdl.sh` doesn't work, the access may be denied by Google Drive because many users downloaded the models within 24 hours. In this case, you can manually download the models with the following links:
* [models1.tar](https://drive.google.com/u/0/uc?id=1az3X7Z88RVFYEAExXK_njJTlDb-CPrc4&export=download) (25.6G, required, the models to predict host taxa from phylum to genus)
* [models2.tar](https://drive.google.com/u/0/uc?id=1ldCUZT5OchCtekUegf-uFABDhdATi-Zg&export=download) (18.3G, optional, the models to predict species)

After the downloading, you can place the file `models1.tar` or the two files under the path `HOTSPOT/models/` and run the Python script `uncompress.py` to uncompress them:
```
python uncompress.py models
```
Then, you have completed the preparation steps and can use HOTSPOT for host prediction.


# Advanced datasets
We collected multiple datasets containing well-annotated plasmid contigs (>1.5kbp) assembled from different metagenomic data to test HOTSPOT's performance in the paper. These datasets and the corresponding taxonomic label information are uploaded to the `Datasets` folder. The labels annotated with `(novel)` refer to the hosts which are rare and do not exist in our label set. Users can use these datasets to become familiar with HOTSPOT or benchmark with other host prediction tools. The datasets comprise:
* [Simulated metagenomic data](https://github.com/fmaguire/MAG_gi_plasmid_analysis),
* Mock metagenomic data (four datasets with accession number SRR072232, SRR072233, SRR172902, and SRR172903),
* [The Hi-C dataset](https://osf.io/ezb8j/wiki/home/),
* [The CAMI2 marine S0 dataset](https://www.microbiome-cosi.org/cami/cami/cami2).

The detailed preprocessing steps of these datasets are described in our paper.


# Usage
Before predicting the hosts, you have to run `preprocessing.py`, which filters lengths and encodes features for input plasmid sequences. Then, you can use `HOTSPOT.py` for host prediction with the pre-trained Transformer models. By default, the temporary files are stored in the folder `temporary_files/`, and the prediction results are stored in the TSV file `Results/host_lineage.tsv`.

## Simple example
```
python preprocessing.py --contigs Example_fasta/multiple_plasmids.fasta
python HOTSPOT.py       # Recommend using gpu to accelerate the program
```


## Running time evaluation
When using the HOTSPOT tool, only one pre-trained classifier is loaded into the GPU at a time, enabling successful running even on a computer with a small GPU memory. In addition, despite the large model size, the prediction process of HOTSPOT is fast. For example, the running time of HOTSPOT on the complete plasmid test set (7,186 plasmids, 530M, for more detail, please see the HOTSPOT paper) with 8 threads is listed in the table below:

| preprocessing.py | HOTSPOT.py | Total running time |
| ------------- | ------------- | ------------- |
| 6.7795h | 3.0805min | 6.8308h |

Notably, the majority of the running time is used to run Prodigal and DIAMOND BLASTP.


## Early stop mechanism using the Monte Carlo dropout (MC-dropout)
HOTSPOT provides two special modes, *specific mode* and *accurate mode*, aiming at higher prediction accuracy using the MC-dropout based early stop of the tree search. To enable the early stop, you can use either the option `--mode 2` (*specific mode*) or `--mode 3` (*accurate mode*) when running `HOTSPOT.py`. Notably, the *accurate mode* has a more stringent uncertainty cutoff than the *specific mode*, leading to more accurate prediction but returning taxa in higher levels for some inputs. In addition, the number of the dropout-enabled forward passes (using variance to estimate the prediction uncertainty) can be chosen by the option `--mcnum` (default: 100).

For example:
```
python HOTSPOT.py --mode 3 --mcnum 80
```


## Full command-line options
preprocessing.py:
```
The usage of preprocessing.py:
            [-h, --help]          Show the help message and exit
            [--contigs INPUT_FA]  FASTA file of the inputs (one or more contigs, default test_contigs.fa)
            [--len MINIMUM_LEN]   Minimum length (bp) of contigs for length filter (default 1500)
            [--threads NUM]       Number of threads to use (default 8)
            [--dbdir DR]          Database directory (default database/)
            [--midfolder DIR]     Folder to store the intermediate files from preprocessing (default temporary_files/)
```

HOTSPOT.py:
```
The usage of HOTSPOT.py:
            [--midfolder DIR]   Folder to store the intermediate files from preprocessing (default temporary_files/)
            [--mdldir DR]       Pre-trained models' directory (default models/)
            [--dbdir DR]        Database directory (default database/)
            [--out OUT]         Path of the output file (default Result/prediction.tsv)
            [--threads NUM]     Number of threads to use if 'cpu' is detected ('cuda' not found, default 8)
            [--mode MOD]        Selected early stop mode.
                                1: sensitive mode (no early stop used, default)
                                2: specific mode (enabling the early stop)
                                3: accurate mode (enabling the early stop with more stringent uncertainty cutoff, leading to more accurate prediction but returning taxa in higher levels for some inputs)
                                (default 1)
            [--mcnum MC]        The number of the dropout-enabled forward passes to estimate the uncertainty (default: 100, minimum: 10)
```


# Format of the output file
1. The format of the input file should be a fasta file that contains one or more plasmid contigs or complete plasmids. If an input contig is too short for protein translation by Prodigal (e.g., <1.5kbp), HOTSPOT may not output any result for this contig.
2. The output is a TSV file containing the predicted host information from phylum to species. Each row corresponds to one input contig. Examples:

| Contig | phylum | class | order | family | genus | species |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| CP090713.1  | Proteobacteria  | Betaproteobacteria  | Burkholderiales  | Burkholderiaceae  | Burkholderia  | Burkholderia multivorans  |
| NZ_CP050042.1  | Proteobacteria  | Gammaproteobacteria  | Enterobacterales  | Enterobacteriaceae  | Escherichia  | Escherichia coli  | 
| NZ_CP083619.1  | Firmicutes  | Clostridia  | Eubacteriales  | Peptostreptococcaceae  | Clostridioides  | Clostridioides difficile  |
| NZ_CP083659.1  | Proteobacteria  | Gammaproteobacteria  | Moraxellales  | Moraxellaceae  | Acinetobacter  | Acinetobacter variabilis  |
| Z22927.1  | Actinobacteria  | Actinomycetia  | Corynebacteriales  | Corynebacteriaceae  | Corynebacterium  | Corynebacterium glutamicum  |


# Contact
If you have any questions, please email us:
  
>yongxinji2-c@my.cityu.edu.hk (Yongxin JI)
  
>jyshang2-c@my.cityu.edu.hk (Jiayu SHANG)
  
>xubotang2-c@my.cityu.edu.hk (Xubo TANG)
