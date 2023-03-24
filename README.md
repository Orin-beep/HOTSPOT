# HOTSPOT: Hierarchical hOst predicTion for aSsembled Plasmid cOntigs with Transformer

HOTSPOT is a learning-based tool designed for plasmid host prediction. Its backbone is a phylogenetic tree of plasmids' hosts (bacteria) from phylum to species. The top-down tree search can accurately predict the hosts' taxonomic labels by incorporating the state-of-the-art language model, Transformer, in each node’s taxon classifier. __To use HOTSPOT, you only need to input plasmid sequences (complete or segmented) into the program.__


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
* For gpu version pytorch: Search [pytorch](https://pytorch.org/get-started/previous-versions/) to find the correct cuda version according to your computer


## Prepare the environment
After cloning this repository (```git clone https://github.com/Orin-beep/HOTSPOT```), you can use anaconda to install the environment.yaml. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Otherwise, it will run with cpu version). The command is: 
```
cd HOTSPOT/
conda env create -f environment.yaml -n hotspot
conda activate hotspot
```


## Prepare the database and pre-trained models (from Google Drive)
To download the database and pre-trained models, you can simply use these bash scripts: 
```
sh prepare_db.sh        # download and unzip the database, 161.2 MB
sh prepare_mdl.sh       # download and unzip the models, 13.14 GB
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


# Test datasets
You can download the plasmid contigs and raw data of the datasets evluated in the paper 'HOTSPOT: Hierarchical hOst predicTion for aSsembled Plasmid cOntigs with Transformer' with the following links:
| Datasets | Annotated plasmid contigs | Raw data and description |
| ------------- | ------------- | ------------- |
| Simulated metagenomic data | [plasmid_contigs_mag.fa](https://drive.google.com/file/d/1jUyKpcJZj2gLJtnuuj-UoD4Wsd1wYdT6/view?usp=share_link) | [original_dataset.tar.gz](https://drive.google.com/file/d/1fZgyb5jHmSwMKzB3Co4mVf_AMFVsurDJ/view) contains the assembled contigs and the code generating the simulated data|
| Mock metagenomic data | [SRR072232.fasta](https://drive.google.com/file/d/1fMcbmGkMfPBfx1T2ysEVL-uCoxwCFbxQ/view), [SRR072233.fasta](https://drive.google.com/file/d/1xmgfTFkbcu0neZS_qVqmqOpS6Lhjh9FG/view), [SRR172902.fasta](https://drive.google.com/file/d/1Ya548WWe9UPXNgBsqfleap0n1u4xbamq/view), [SRR172903.fasta](https://drive.google.com/file/d/1l6msS2l7qGjYOWAqlMNUE7vOw5w2Hjbm/view) | [SRR072232](https://www.ncbi.nlm.nih.gov/sra/SRR072232), [SRR072233](https://www.ncbi.nlm.nih.gov/sra/?term=SRR072233), [SRR172902](https://www.ncbi.nlm.nih.gov/sra/?term=SRR172902), [SRR172903](https://www.ncbi.nlm.nih.gov/sra/?term=SRR172903). The reference genomes: [reference_genomes_mock.fasta](https://drive.google.com/file/d/159VNsSoPbMQWJmQnPYqqW4lJ8bTLTN6l/view) |
| Hi-C dataset | [plasmid_contigs_hi-c.fa](https://drive.google.com/file/d/1ZysHSuT2afdqDCkYxs44XCk4aqyYXFyu/view?usp=share_link) | [wastewater_hi-c_data.tar.gz](https://drive.google.com/file/d/1Ka_eMr9LDkNU5nZJW7OP6Ea1OeiWFuM2/view). Data source: https://osf.io/ezb8j/wiki/home/ |


# Running time evaluation
We run HOTSPOT with 8 threads and gpu on 4,536 complete plasmids (333MB) as an example. The required running time for the two steps is listed below:

| preprocessing.py | HOTSPOT.py | Total running time |
| ------------- | ------------- | ------------- |
| 3h38min | 52.65s | 3h39min |

Thus, most of the time is used to run Prodigal and DIAMOND BLASTP for preprocessing.


# Contact
If you have any questions, please email us:
  
>yongxinji2-c@my.cityu.edu.hk (Yongxin JI)
  
>jyshang2-c@my.cityu.edu.hk (Jiayu SHANG)
  
>xubotang2-c@my.cityu.edu.hk (Xubo TANG)
