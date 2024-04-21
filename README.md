# HOTSPOT

HOTSPOT is a deep learning tool based on the Transformer model, specifically designed for predicting the hosts of plasmids. __To use HOTSPOT, you can easily input your plasmid DNA sequences (complete plasmids or contigs) into the program.__ Subsequently, you will obtain prediction results of the hosts' taxa, ranging from phylum to species level.


### E-mail: yongxinji2-c@my.cityu.edu.hk


#### *__[Update - 2024 - 04 - 21]__* :  <BR/>
* *We substantially reduce HOTSPOT's model size by training six models, with each model corresponding to one of the six taxonomic levels. <BR/>*


# Install (Linux or ubuntu only)
## Dependencies
* [Python 3.x](https://www.python.org/downloads/)
* [NumPy](https://pypi.org/project/numpy/)
* [bidict](https://pypi.org/project/bidict/)
* [Pandas](https://pypi.org/project/pandas/)
* [PyTorch](https://pytorch.org/get-started/previous-versions/)>1.8.0
* [diamond](https://anaconda.org/bioconda/diamond)
* [Prodigal](https://anaconda.org/bioconda/prodigal)
* [biopython](https://pypi.org/project/biopython/)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)=2.13.0 (conda install -c bioconda blast=2.13.0)
* [HMMER](https://anaconda.org/bioconda/hmmer)
* [treelib](https://pypi.org/project/treelib/)

If you want to use the GPU to accelerate the program:
* CUDA
* PyTorch-GPU
* For CPU version PyTorch: ```conda install pytorch torchvision torchaudio cpuonly -c pytorch```
* For GPU version PyTorch: search [pytorch](https://pytorch.org/get-started/previous-versions/) to find the correct CUDA version according to your computer


## Prepare the environment
After cloning this repository (```git clone https://github.com/Orin-beep/HOTSPOT```), you can use Anaconda to install ```environment.yaml```. This will install all packages you need in GPU mode (make sure you have installed CUDA on your system to use the GPU version; otherwise, HOTSPOT will run in CPU mode). The installing command is: 
```
cd HOTSPOT/
conda env create -f environment.yaml -n hotspot
conda activate hotspot
```


## Prepare default database and models (from Google Drive)
To download the default database and models, you can use the following bash scripts: 
```
sh prepare_db.sh        # download and unzip database, 161.2 MB
sh prepare_mdl.sh       # download and unzip models, 13.14 GB
```


## Download the database and pre-trained models manually
If the bash scripts do not work, you can manually download the database and models using the following links:
* [database.tgz](https://drive.google.com/file/d/1VIbfp35X5JMiA7BfOS3lNBycT73uFcTN/view)
* [models.tgz](https://drive.google.com/file/d/1L6ogZhdAWJ7Ns8Hz59m7W2kFPMcMTC6u/view)

After downloading the `database.tgz` and `models.tgz` to HOTSPOT's main directory, you have to unzip them:
```
tar -zxvf database.tgz
rm database.tgz

tar -zxvf models.tgz
rm models.tgz
```


# Usage
Before predicting the hosts, you have to run `preprocessing.py`, which filters lengths and encodes features for input plasmid sequences. Then, you can use `HOTSPOT.py` for host prediction with the pre-trained Transformer models. By default, the temporary files are stored in the folder `temporary_files/`, and the prediction results are stored in the TSV file `Results/host_lineage.tsv`.

## Simple example
```
python preprocessing.py --contigs Example_fasta/multiple_plasmids.fasta
python HOTSPOT.py       # Recommend using gpu to accelerate the program
```


# Format of the output file
The output is a TSV file containing the predicted host lineages from phylum to species. Each row corresponds to one input plasmid contig. For example:

| Contig | phylum | class | order | family | genus | species |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| NZ_CP050042.1  | Pseudomonadota  | Gammaproteobacteria  | Enterobacterales  | Enterobacteriaceae  | Escherichia  | Escherichia coli  | 
| NZ_CP083619.1  | Bacillota  | Clostridia  | Eubacteriales*  | Peptostreptococcaceae  | Clostridioides*  | Clostridioides difficile*  |
| NZ_CP083659.1  | Pseudomonadota  | Gammaproteobacteria  | Moraxellales  | Moraxellaceae*  | Acinetobacter  | Acinetobacter variabilis  |
| Z22927.1  | Actinomycetota  | Actinomycetes*  | Corynebacteriales  | Corynebacteriaceae  | Corynebacterium*  | Corynebacterium glutamicum*  |

Notably, the taxon labeled with a star `*` is not predicted by the taxon classifier because its parent node has only one child in the tree. 
The current phylogenetic tree used by HOTSPOT is smaller than the complete bacterial phylogenetic tree because: 1) not all bacteria contain plasmids, and 2) the host taxa covered by available sequenced plasmids are limited. Thus, we advise users to examine starred taxa more carefully.


## Full command-line options
preprocessing.py:
```
The usage of preprocessing.py:
            [-h, --help]          Show the help message and exit
            [--contigs INPUT_FA]  FASTA file of the input sequences (one or more contigs in a single FASTA file, default test_contigs.fa)
            [--len MINIMUM_LEN]   Minimum length (bp) of contigs for length filtering (default 1500)
            [--threads NUM]       Number of threads to use (default 8)
            [--dbdir DR]          Database directory (default database/)
            [--midfolder DIR]     Folder to store the intermediate files from preprocessing (default temporary_files/)

```

HOTSPOT.py:
```
The usage of HOTSPOT.py:
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

```

## Early stop mechanism using the Monte Carlo dropout (MC-dropout)
HOTSPOT provides two special modes, *specific mode* and *accurate mode*, aiming at higher accuracy using the MC-dropout based early stop for the tree search. To enable the early stop, you can use the option `--accurate 1` when running `HOTSPOT.py`, and the results of the two modes will be stored in the output directory. Specifically, the *accurate mode* has a more stringent uncertainty cutoff than the *specific mode*, leading to more accurate prediction but returning taxa in higher levels for some inputs. In addition, the number of dropout-enabled forward passes can be chosen by the option `--mcnum` (default: 100).

For example (the prediction will take more time):
```
python HOTSPOT.py --accurate 1
```


