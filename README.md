# HOTSPOT: plasmid host prediction

HOTSPOT's backbone is a phylogenetic tree of the plasmid hosts (bacteria) from phylum to species. By incorporating the state-of-the-art language model, Transformer, in each node’s taxon classifier, the top-down tree search can accurately predict the host taxonomy for the input plasmid contigs. There are totally 115 taxon classifiers, each corresponding to a node with more than one child node. __To use HOTSPOT, you only need to input complete plasmids or plasmid contigs assembled from metagenomic data into the program.__


# Required Dependencies
* Python 3.x
* Numpy
* Pandas
* Pytorch>1.8.0
* [Diamond](https://github.com/bbuchfink/diamond)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Biopython](https://biopython.org/)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [HMMER](http://hmmer.org/)
* [treelib](https://pypi.org/project/treelib/)
* [gdown](https://github.com/wkentaro/gdown) (to download the database and pre-trained models from Google Drive)

If you want to use the gpu to accelerate the program:
* cuda
* Pytorch-gpu
* For cpu version pytorch: ```conda install pytorch torchvision torchaudio cpuonly -c pytorch```
* For gpu version pytorch: Search [pytorch](https://pytorch.org/) to find the correct cuda version according to your computer


## Prepare the environment
After cloning this respository, you can use anaconda to install the environment.yaml. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: 
```
conda env create -f environment.yaml -n hotspot
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

After the downloading, you can place the file `models1.tar` or the two files under the path `HOTSPOT/models` and run the Python script `uncompress.py` to uncompress them:
```
python uncompresse.py models
```
Then, you have completed the preparation steps and can use HOTSPOT for host prediction.


# Usage
Before running the host prediction, you have to use __`preprocessing.py` for length filter and feature encoding__ (including the hypothetical proteins, identified MOB/MPF proteins and Inc group of the inputs. Details are available in the paper "HOTSPOT: Hierarchical hOst predicTion for aSsembled Plasmid cOntigs with Transformer"). The temporary files will be stored in the folder `temporary_files/` by default. Then, `HOTSPOT.py` performs the host prediction with the phylogenetic tree search and pre-trained Transformer models. The results are stored in the TSV file `Result/prediction.tsv` by default.


## Simple example
```
python preprocessing.py --contigs Example_fasta/multiple_plasmids.fasta
python HOTSPOT.py       # It would highly recommend using gpu to accelerate the program
```


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
