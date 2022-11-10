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


## Prepare the database
To prepare the database, you can simply use this bash script: 
```
sh prepare_db.sh
```
To be sure that the downloaded `database/` folder is placed in the HOTSPOT's main directory. Otherwise, you have to use option `--dbdir` to specify the database path when running `preprocessing.py` and `HOTSPOT.py`.


## Prepare the pre-trained models
If you don't need the host identification at the species level, you may only download the pre-trained models from phylum to genus (26G):
```
sh prepare_light_mdl.sh
```

Or if you have sufficient disk storage, you can download the complete models (45G):
```
sh prepare_mdl.sh
```

Similar to the database preparation, if the model path is not under HOTSPOT's main directory, you need to use option `--mdldir` to specify the model path when running `HOTSPOT.py`.


The pre-trained models are a little big because the phylogenetic includes 10, 19, 46, 94, 180, and 272 nodes for the six taxonomic ranks from phylum to species, respectively. We will try to optimize the model size in the next version of HOTSPOT.

# Usage
Before the prediction, you have to use `preprocessing.py` for length filter and sentence matrix generation (as the input into the Transformer model). The temporary files will be stored in the folder `temporary_files/` by default. Then, `HOTSPOT.py` will help you with the host prediction. The results are stored in the csv file `Result/prediction.tsv` by default.

## Simple example
```
python preprocessing.py --contigs Example_fasta/multiple_plasmids.fasta
python HOTSPOT.py
```

## Early stop function using Monte Carlo (MC) dropout
HOTSPOT provides two special modes, `specific mode` and `accurate mode`, aiming to achieve higher precision and handle broad-host-range (BHR) plasmids. To enable these two modes, you can use options `--mode 2` or `--mode 3`. In addition, the prediction number for MC dropout can also be chosen by the option `--mcnum`.

For example:
```
python HOTSPOT.py --mode 3 --mcnum 80
```


## Full command-line options
preprocessing.py:
```
The usage of preprocessing.py:
            [--contigs INPUT_FA]  Input fasta file (containing one or multiple plasmid contigs)
            [--len MINIMUM_LEN]   Predict only for sequence >= len bp (default 1500)
            [--threads NUM]       Number of threads to run preprocessing (default 8)
            [--dbdir DR]          Path to store the database directory (default database/)
            [--midfolder DIR]     Folder to store the intermediate files (default temporary_files/)

```
HOTSPOT.py:
```
The usage of HOTSPOT.py:
            [--midfolder DIR]   Intermediate file folder output by preprocessing.py (default temporary_files/)
            [--mdldir DR]       Path to store the HOTSPOT pre-trained model directory (default models/)
            [--dbdir DR]        Path to store the database directory (default database/)
            [--out OUT]         File path to store the prediction results (default "Result/prediction.tsv")
            [--threads NUM]     Number of threads to run if 'cpu' is detected ('cuda' not found) (default 8)
            [--mode MOD]        Three modes with different Monte Carlo (MC) dropout early stop parameters. If 2 or 3 is chosen, the prediction process will slow down.
                                1: sensitive mode (lower precision, higher recall)  (default)
                                2: specific mode (moderate precision, moderate recall)
                                3: accurate mode (higher precision, lower recall)
                                (default 1)
            [--mcnum MC]        The number of MC dropout predictions (The minimum value is 10). If you enable the early stop mode (2 or 3), the running time will increase in proportion to this number (default 100)
```


# Format of the output file
1. The format of the input file should be a fasta file that contains one or more plasmid contigs (complete plasmids are also ok). If an input contig is too short (e.g., <1.5kbp), it's possible that **Prodigal** can predict no protein. In this case, HOTSPOT will not output any result for this contig.
2. The output file format is a csv file containing the predicted host taxonomic lineage information of input plasmid contigs. Each row represents the prediction result for one input contig. **Contig** is the contig accession number from the input.

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
