# HOTSPOT

HOTSPOT is a deep learning tool based on the Transformer model, specifically designed for predicting the hosts of plasmids. __To use HOTSPOT, you can easily input your plasmid DNA sequences (complete plasmids or contigs) into the program.__ Subsequently, you will obtain prediction results of the hosts' taxa, ranging from phylum to species level.


### E-mail: yongxinji2-c@my.cityu.edu.hk


### Version: V1.1 (update at 2024-04-21)

#### *__[Update - 2024 - 04 - 21]__* :  <BR/>
* *We substantially reduce HOTSPOT's model size by training six models, with each model corresponding to one of the six taxonomic levels. Furthermore, we have incorporated an additional feature that allows you to train your custom models using your own plasmid dataset. <BR/>*


# Install (Linux or Ubuntu only)
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

If you want to use GPU to accelerate the program:
* CUDA
* PyTorch-GPU
* For CPU version PyTorch: ```conda install pytorch torchvision torchaudio cpuonly -c pytorch```
* For GPU version PyTorch: search [PyTorch](https://pytorch.org/get-started/previous-versions/) to find the correct CUDA version according to your computer


## Prepare the environment
After cloning this repository (```git clone https://github.com/Orin-beep/HOTSPOT```), you can use Anaconda to install ```environment.yaml```. This will install all packages you need in GPU mode (make sure you have installed CUDA on your system to use the GPU version; otherwise, HOTSPOT will run in CPU mode). The installing command is: 
```
git clone https://github.com/Orin-beep/HOTSPOT
cd HOTSPOT/
conda env create -f environment.yaml -n hotspot
conda activate hotspot
```
If Anaconda fails to work, you can prepare the environment by individually installing the packages listed in the __Dependencies__ section.

## Prepare default database and models (from Google Drive)
To download the default database and models, you can use the following bash scripts: 
```
sh prepare_db.sh        # download and unzip the database folder, 91.8 MB
sh prepare_mdl.sh       # download and unzip the model folder, 1.78 GB
```


## Alternative way: download default database and models manually
If the above bash scripts do not work, you can manually download the default database and models using the following Google Drive links:
* [database.tar.gz](https://drive.google.com/file/d/1ZSTz3kotwF8Zugz_aBGDtmly8BVo9G4T/view)
* [models.tar.gz](https://drive.google.com/file/d/1bnA1osvYDgYBi-DRFkP-HrvcnvBvbipF/view)

After downloading the `database.tar.gz` and `models.tar.gz` to HOTSPOT's main directory, you have to unzip them:
```
tar -zxvf database.tar.gz
rm database.tar.gz

tar -zxvf models.tar.gz
rm models.tar.gz
```


# Usage
Before running HOTSPOT, you should run `preprocessing.py` to encode the input plasmid sequences into sentences. After that, you can use `hotspot.py` for host prediction.

## Simple example
```
python preprocessing.py --fasta example_plasmids/NZ_CP083659.fasta --database database --model_path models
python hotspot.py
```


## Format of the output file
The results will be saved in a TSV file (default: results/host_lineage.tsv) containing the predicted host lineages from phylum to species level. Each row corresponds to an input plasmid sequence. Examples:

| Contig | Phylum | Class | Order | Family | Genus | Species |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| NZ_CP050042.1  | Pseudomonadota  | Gammaproteobacteria  | Enterobacterales  | Enterobacteriaceae  | Escherichia  | -  | 
| NZ_CP083619.1  | Bacillota  | Clostridia  | Eubacteriales  | Peptostreptococcaceae  | Clostridioides  | Clostridioides difficile  |
| NZ_CP083659.1  | Pseudomonadota  | Gammaproteobacteria  | Moraxellales  | Moraxellaceae  | Acinetobacter  | Acinetobacter variabilis  |
| Z22927.1  | Actinomycetota  | Actinomycetes  | Corynebacteriales  | Corynebacteriaceae  | Corynebacterium  | Corynebacterium glutamicum  |

The dash '-' indicates that HOTSPOT cannot provide accurate predictions for the corresponding input at this taxonomic level.


## Accurate mode with Monte Carlo dropout (MC-dropout)
HOTSPOT provides a specialized *accurate mode* designed to enhance accuracy using MC-dropout based early stop. You can activate this mode using the `--accurate_mode True` option. For advanced users, you can utilize the `--mc_num`, `--min_mean`, and `--max_var` options to have control over the specific implementation details of the MC-dropout mechanism. Detailed usage of these options is shown in the following __Full command-line options__ section. Example:

```
python hotspot.py --accurate_mode True
```


## Train your custom models



## Full command-line options
preprocessing.py:
```
Usage of preprocessing.py:
        [--fasta FASTA] FASTA file of the plasmid DNA sequences to be predicted (either complete sequences or contigs), default: multiple_plasmids.fasta
        [--database DATABASE]   path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: database
        [--model_path MODEL_PATH]   path of the folder storing the downloaded or your customized models, default: models
        [--midfolder]   folder to store the intermediate files for prediction, default: temp
        [--len LEN] minimum length of plasmid DNA sequences, default: 1500
        [--threads THREADS] number of threads utilized for preprocessing, default: 2
```

hotspot.py:
```
Usage of hotspot.py:
        [--model_path MODEL_PATH]   path of the folder storing the downloaded or your customized models, default: models
        [--midfolder MIDFOLDER] folder to store the intermediate files (generated by preprocessing.py), default: temp
        [--device DEVICE]   device utilized for prediction ('gpu' or 'cpu'), default: 'gpu'
        [--threads THREADS] number of threads utilized for prediction if 'cpu' is detected ('cuda' not found), default: 2
        [--batch_size BATCH_SIZE]   batch size for prediction, default: 200
        [--out OUT] path to store the prediction results, default: results
        [--accurate_mode ACCURATE_MODE] if this option is set to True, HOTSPOT will run in accurate mode. In accurate mode, the prediction process is slightly slower due to the activation of the MC-dropout mechanism. Additionally, besides the normal output file 'host_lineage.tsv', a supplementary TSV file named 'host_lineage_acc.tsv' will be generated. This file contains high-confidence predictions with a sacrifice of resolution, meaning that for certain inputs, the returned taxa may be at higher levels in the taxonomic hierarchy, default: False
        [--mc_num MC_NUM]   if the accurate mode is activated, you can use this option to specify the number of dropout-enabled forward passes for the MC-dropout mechanism, default: 100
        [--min_mean MIN_MEAN]   the minimum mean value for a prediction that will not trigger early stopping by MC-dropout, default: 0.75
        [--max_var MAX_VAR] the maximum variance value for a prediction that will not trigger early stopping by MC-dropout, default: 0.3
```

preprocessing_train.py:
```
Usage of preprocessing_train.py:
        [--fasta FASTA]    FASTA file of plasmid DNA sequences for training (preferably complete sequences, default: '../training_dataset/plasmids.fasta')
        [--host_info HOST_INFO] TSV file containing complete host lineage information from phylum to species for each training plasmid, default: '../training_dataset/host_lineages.tsv'
        [--database DATABASE]    path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: '../database'
        [--model_path MODEL_PATH]  folder to store your customized models, default: 'models'
        [--midfolder MIDFOLDER]   folder to store the intermediate files, default: 'temp'
        [--len LEN] minimum length of plasmid DNA sequences, default: 1500
        [--train_val_list TRAIN_VAL_LIST]  TXT file containing the information of the training/validation sets. The first row should display the list of training plasmids, while the second row should display the list of validation plasmids. Each plasmid in the lists should be separated by a space, default: '../training_dataset/train_val.txt'
        [--train_ratio TRAIN_RATIO] the ratio of the training set size to the total number of input training plasmids. If the train_val_list file is not provided, the training plasmids will be randomly split into train/validation sets using the specified ratio, default: 0.8
        [--labels LABELS]  TXT file containing the specified taxonomic labels for different levels. The file should comprise six rows, where each row corresponds to the labels for the phylum, class, order, family, genus, and species levels, respectively. Within each row, the labels should be separated by tabs ('	'), default: None
        [--num_plasmids NUM_PLASMIDS]    minimun number of training plasmids associated with a taxonomic label (used when the labels file is not provided), default: 20
        [--add_frags ADD_FRAGS]   whether to augment the training set by randomly cutting fragments ranging from 1.5 to 15 kbp (this may slightly slow down the training process, but it will significantly enhance the performance of host prediction), default: False
        [--num_frags NUM_FRAGS]   maximum number of added fragments from each training plasmid, default: 5
        [--threads THREADS] number of threads utilized for preprocessing, default: 2
```

train.py
```
Usage of train.py:
        [--model_path MODEL_PATH]   folder to store your customized models, default: models
        [--midfolder MIDFOLDER] folder to store the intermediate files (generated by preprocessing_train.py), default: temp
        [--device DEVICE]   device utilized for training ('gpu' or 'cpu'), default: 'gpu'
        [--threads THREADS] number of threads utilized for training if 'cpu' is detected ('cuda' not found), default: 2
        [--lr LR]   learning rate for training the models, default: 0.005
        [--batch_size BATCH_SIZE]   batch size for training the models, default: 200
        [--epoch_num EPOCH_NUM] number of epochs for training the models, default: 20
        [--dropout DROPOUT] dropout rate for training the models, default: 0.5 
```


## References
How to cite this tool:
```
Ji, Yongxin, et al. "HOTSPOT: hierarchical host prediction for assembled plasmid contigs with transformer." Bioinformatics 39.5 (2023): btad283. https://academic.oup.com/bioinformatics/article/39/5/btad283/7136643
```


## The default plasmid hosts' phylogenetic tree
```
Bacteria
├── Actinomycetota
│   └── Actinomycetes
│       ├── Bifidobacteriales
│       │   └── Bifidobacteriaceae
│       │       └── Bifidobacterium
│       │           └── Bifidobacterium longum
│       ├── Kitasatosporales
│       │   └── Streptomycetaceae
│       │       └── Streptomyces
│       ├── Micrococcales
│       │   ├── Microbacteriaceae
│       │   │   ├── Clavibacter
│       │   │   │   └── Clavibacter michiganensis
│       │   │   └── Curtobacterium
│       │   └── Micrococcaceae
│       │       ├── Arthrobacter
│       │       └── Paenarthrobacter
│       ├── Mycobacteriales
│       │   ├── Corynebacteriaceae
│       │   │   └── Corynebacterium
│       │   │       └── Corynebacterium glutamicum
│       │   ├── Gordoniaceae
│       │   │   └── Gordonia
│       │   ├── Mycobacteriaceae
│       │   │   ├── Mycobacterium
│       │   │   │   ├── Mycobacterium avium
│       │   │   │   └── Mycobacterium intracellulare
│       │   │   ├── Mycobacteroides
│       │   │   │   └── Mycobacteroides abscessus
│       │   │   └── Mycolicibacterium
│       │   ├── Nocardiaceae
│       │   │   ├── Nocardia
│       │   │   ├── Prescottella
│       │   │   │   └── Prescottella equi
│       │   │   └── Rhodococcus
│       │   │       ├── Rhodococcus erythropolis
│       │   │       ├── Rhodococcus pyridinivorans
│       │   │       └── Rhodococcus qingshengii
│       │   └── Tsukamurellaceae
│       │       └── Tsukamurella
│       │           └── Tsukamurella tyrosinosolvens
│       ├── Propionibacteriales
│       ├── Pseudonocardiales
│       │   └── Pseudonocardiaceae
│       └── Streptosporangiales
├── Bacillota
│   ├── Bacilli
│   │   ├── Bacillales
│   │   │   ├── Bacillaceae
│   │   │   │   ├── Bacillus
│   │   │   │   │   ├── Bacillus anthracis
│   │   │   │   │   ├── Bacillus cereus
│   │   │   │   │   ├── Bacillus cytotoxicus
│   │   │   │   │   ├── Bacillus licheniformis
│   │   │   │   │   ├── Bacillus mycoides
│   │   │   │   │   ├── Bacillus paranthracis
│   │   │   │   │   ├── Bacillus subtilis
│   │   │   │   │   ├── Bacillus thuringiensis
│   │   │   │   │   ├── Bacillus toyonensis
│   │   │   │   │   └── Bacillus velezensis
│   │   │   │   ├── Geobacillus
│   │   │   │   └── Priestia
│   │   │   │       ├── Priestia aryabhattai
│   │   │   │       └── Priestia megaterium
│   │   │   ├── Listeriaceae
│   │   │   │   └── Listeria
│   │   │   │       ├── Listeria innocua
│   │   │   │       └── Listeria monocytogenes
│   │   │   ├── Paenibacillaceae
│   │   │   │   └── Paenibacillus
│   │   │   ├── Planococcaceae
│   │   │   │   └── Planococcus
│   │   │   └── Staphylococcaceae
│   │   │       ├── Macrococcus
│   │   │       │   └── Macrococcus caseolyticus
│   │   │       ├── Mammaliicoccus
│   │   │       │   └── Mammaliicoccus sciuri
│   │   │       └── Staphylococcus
│   │   │           ├── Staphylococcus aureus
│   │   │           ├── Staphylococcus capitis
│   │   │           ├── Staphylococcus epidermidis
│   │   │           ├── Staphylococcus equorum
│   │   │           ├── Staphylococcus haemolyticus
│   │   │           ├── Staphylococcus hominis
│   │   │           ├── Staphylococcus lugdunensis
│   │   │           ├── Staphylococcus saprophyticus
│   │   │           ├── Staphylococcus simulans
│   │   │           └── Staphylococcus warneri
│   │   └── Lactobacillales
│   │       ├── Carnobacteriaceae
│   │       │   └── Carnobacterium
│   │       ├── Enterococcaceae
│   │       │   ├── Enterococcus
│   │       │   │   ├── Enterococcus casseliflavus
│   │       │   │   ├── Enterococcus faecalis
│   │       │   │   ├── Enterococcus faecium
│   │       │   │   └── Enterococcus hirae
│   │       │   └── Vagococcus
│   │       ├── Lactobacillaceae
│   │       │   ├── Companilactobacillus
│   │       │   ├── Lacticaseibacillus
│   │       │   │   ├── Lacticaseibacillus paracasei
│   │       │   │   └── Lacticaseibacillus rhamnosus
│   │       │   ├── Lactiplantibacillus
│   │       │   │   ├── Lactiplantibacillus pentosus
│   │       │   │   └── Lactiplantibacillus plantarum
│   │       │   ├── Lactobacillus
│   │       │   ├── Latilactobacillus
│   │       │   │   ├── Latilactobacillus curvatus
│   │       │   │   └── Latilactobacillus sakei
│   │       │   ├── Lentilactobacillus
│   │       │   ├── Leuconostoc
│   │       │   │   ├── Leuconostoc citreum
│   │       │   │   └── Leuconostoc mesenteroides
│   │       │   ├── Levilactobacillus
│   │       │   │   └── Levilactobacillus brevis
│   │       │   ├── Ligilactobacillus
│   │       │   │   └── Ligilactobacillus salivarius
│   │       │   ├── Limosilactobacillus
│   │       │   │   ├── Limosilactobacillus fermentum
│   │       │   │   └── Limosilactobacillus reuteri
│   │       │   ├── Loigolactobacillus
│   │       │   │   └── Loigolactobacillus backii
│   │       │   ├── Pediococcus
│   │       │   │   ├── Pediococcus damnosus
│   │       │   │   └── Pediococcus pentosaceus
│   │       │   ├── Secundilactobacillus
│   │       │   └── Weissella
│   │       │       └── Weissella cibaria
│   │       └── Streptococcaceae
│   │           ├── Lactococcus
│   │           │   ├── Lactococcus cremoris
│   │           │   ├── Lactococcus garvieae
│   │           │   └── Lactococcus lactis
│   │           └── Streptococcus
│   │               ├── Streptococcus suis
│   │               └── Streptococcus thermophilus
│   └── Clostridia
│       └── Eubacteriales
│           ├── Clostridiaceae
│           │   └── Clostridium
│           │       ├── Clostridium botulinum
│           │       ├── Clostridium estertheticum
│           │       └── Clostridium perfringens
│           └── Peptostreptococcaceae
│               └── Clostridioides
│                   └── Clostridioides difficile
├── Bacteroidota
│   ├── Bacteroidia
│   │   └── Bacteroidales
│   │       ├── Bacteroidaceae
│   │       │   ├── Bacteroides
│   │       │   │   ├── Bacteroides fragilis
│   │       │   │   └── Bacteroides thetaiotaomicron
│   │       │   └── Phocaeicola
│   │       └── Tannerellaceae
│   │           └── Parabacteroides
│   │               └── Parabacteroides distasonis
│   ├── Cytophagia
│   │   └── Cytophagales
│   │       ├── Flammeovirgaceae
│   │       │   └── Chondrinema
│   │       │       └── Chondrinema litorale
│   │       ├── Hymenobacteraceae
│   │       │   └── Hymenobacter
│   │       └── Spirosomaceae
│   └── Flavobacteriia
│       └── Flavobacteriales
│           ├── Flavobacteriaceae
│           │   └── Flavobacterium
│           └── Weeksellaceae
├── Campylobacterota
│   └── Epsilonproteobacteria
│       └── Campylobacterales
│           ├── Arcobacteraceae
│           ├── Campylobacteraceae
│           │   └── Campylobacter
│           │       ├── Campylobacter coli
│           │       ├── Campylobacter fetus
│           │       └── Campylobacter jejuni
│           └── Helicobacteraceae
│               └── Helicobacter
│                   └── Helicobacter pylori
├── Chlamydiota
│   └── Chlamydiia
│       └── Chlamydiales
│           └── Chlamydiaceae
│               └── Chlamydia
├── Cyanobacteriota
│   └── Cyanophyceae
│       ├── Chroococcales
│       │   ├── Aphanothecaceae
│       │   └── Chroococcaceae
│       ├── Nostocales
│       │   ├── Calotrichaceae
│       │   │   └── Calothrix
│       │   ├── Nostocaceae
│       │   │   └── Nostoc
│       │   └── Tolypothrichaceae
│       │       └── Tolypothrix
│       ├── Oscillatoriales
│       │   └── Microcoleaceae
│       │       └── Planktothrix
│       ├── Pseudanabaenales
│       └── Synechococcales
│           ├── Acaryochloridaceae
│           │   └── Acaryochloris
│           ├── Merismopediaceae
│           │   └── Synechocystis
│           └── Synechococcaceae
│               └── Synechococcus
├── Deinococcota
│   └── Deinococci
│       ├── Deinococcales
│       │   └── Deinococcaceae
│       │       └── Deinococcus
│       └── Thermales
│           └── Thermaceae
│               └── Thermus
│                   └── Thermus thermophilus
├── Fusobacteriota
│   └── Fusobacteriia
│       └── Fusobacteriales
│           ├── Fusobacteriaceae
│           │   └── Fusobacterium
│           └── Leptotrichiaceae
├── Mycoplasmatota
│   └── Mollicutes
│       ├── Acholeplasmatales
│       │   └── Acholeplasmataceae
│       │       └── Candidatus Phytoplasma
│       ├── Entomoplasmatales
│       │   └── Spiroplasmataceae
│       │       └── Spiroplasma
│       │           └── Spiroplasma citri
│       └── Mycoplasmatales
│           └── Mycoplasmataceae
│               └── Mycoplasmopsis
├── Planctomycetota
├── Pseudomonadota
│   ├── Acidithiobacillia
│   │   └── Acidithiobacillales
│   │       └── Acidithiobacillaceae
│   │           └── Acidithiobacillus
│   ├── Alphaproteobacteria
│   │   ├── Hyphomicrobiales
│   │   │   ├── Aurantimonadaceae
│   │   │   ├── Bartonellaceae
│   │   │   │   └── Bartonella
│   │   │   ├── Brucellaceae
│   │   │   │   └── Brucella
│   │   │   ├── Methylobacteriaceae
│   │   │   │   └── Methylobacterium
│   │   │   ├── Methylocystaceae
│   │   │   │   └── Methylocystis
│   │   │   ├── Nitrobacteraceae
│   │   │   │   └── Bradyrhizobium
│   │   │   ├── Phyllobacteriaceae
│   │   │   │   ├── Aminobacter
│   │   │   │   └── Mesorhizobium
│   │   │   ├── Rhizobiaceae
│   │   │   │   ├── Agrobacterium
│   │   │   │   │   ├── Agrobacterium fabacearum
│   │   │   │   │   ├── Agrobacterium fabrum
│   │   │   │   │   ├── Agrobacterium rhizogenes
│   │   │   │   │   └── Agrobacterium tumefaciens
│   │   │   │   ├── Rhizobium
│   │   │   │   │   └── Rhizobium leguminosarum
│   │   │   │   ├── Shinella
│   │   │   │   └── Sinorhizobium
│   │   │   │       └── Sinorhizobium meliloti
│   │   │   └── Stappiaceae
│   │   │       └── Roseibium
│   │   ├── Rhodobacterales
│   │   │   ├── Paracoccaceae
│   │   │   │   ├── Cereibacter
│   │   │   │   │   └── Cereibacter sphaeroides
│   │   │   │   └── Paracoccus
│   │   │   │       ├── Paracoccus marcusii
│   │   │   │       └── Paracoccus yeei
│   │   │   └── Roseobacteraceae
│   │   │       ├── Leisingera
│   │   │       │   ├── Leisingera aquaemixtae
│   │   │       │   └── Leisingera caerulea
│   │   │       ├── Phaeobacter
│   │   │       │   ├── Phaeobacter gallaeciensis
│   │   │       │   ├── Phaeobacter inhibens
│   │   │       │   └── Phaeobacter piscinae
│   │   │       ├── Pseudosulfitobacter
│   │   │       │   └── Pseudosulfitobacter pseudonitzschiae
│   │   │       ├── Roseivivax
│   │   │       ├── Roseovarius
│   │   │       ├── Salipiger
│   │   │       └── Sulfitobacter
│   │   ├── Rhodospirillales
│   │   │   ├── Acetobacteraceae
│   │   │   │   ├── Acetobacter
│   │   │   │   │   └── Acetobacter pasteurianus
│   │   │   │   ├── Acidiphilium
│   │   │   │   ├── Komagataeibacter
│   │   │   │   └── Roseomonas
│   │   │   │       └── Roseomonas mucosa
│   │   │   └── Azospirillaceae
│   │   │       └── Azospirillum
│   │   ├── Rickettsiales
│   │   │   └── Rickettsiaceae
│   │   │       └── Rickettsia
│   │   └── Sphingomonadales
│   │       ├── Sphingomonadaceae
│   │       │   ├── Novosphingobium
│   │       │   ├── Sphingobium
│   │       │   │   └── Sphingobium yanoikuyae
│   │       │   ├── Sphingomonas
│   │       │   └── Sphingopyxis
│   │       └── Zymomonadaceae
│   │           └── Zymomonas
│   │               └── Zymomonas mobilis
│   ├── Betaproteobacteria
│   │   ├── Burkholderiales
│   │   │   ├── Alcaligenaceae
│   │   │   ├── Burkholderiaceae
│   │   │   │   ├── Burkholderia
│   │   │   │   │   ├── Burkholderia contaminans
│   │   │   │   │   ├── Burkholderia gladioli
│   │   │   │   │   └── Burkholderia glumae
│   │   │   │   ├── Cupriavidus
│   │   │   │   ├── Paraburkholderia
│   │   │   │   └── Ralstonia
│   │   │   └── Comamonadaceae
│   │   │       ├── Comamonas
│   │   │       └── Polaromonas
│   │   ├── Neisseriales
│   │   │   └── Neisseriaceae
│   │   │       └── Neisseria
│   │   │           └── Neisseria gonorrhoeae
│   │   └── Nitrosomonadales
│   └── Gammaproteobacteria
│       ├── Aeromonadales
│       │   └── Aeromonadaceae
│       │       └── Aeromonas
│       │           ├── Aeromonas caviae
│       │           ├── Aeromonas dhakensis
│       │           ├── Aeromonas hydrophila
│       │           ├── Aeromonas salmonicida
│       │           └── Aeromonas veronii
│       ├── Alteromonadales
│       │   ├── Pseudoalteromonadaceae
│       │   │   └── Pseudoalteromonas
│       │   └── Shewanellaceae
│       │       └── Shewanella
│       │           └── Shewanella baltica
│       ├── Enterobacterales
│       │   ├── Enterobacteriaceae
│       │   │   ├── Citrobacter
│       │   │   │   ├── Citrobacter braakii
│       │   │   │   ├── Citrobacter freundii
│       │   │   │   └── Citrobacter portucalensis
│       │   │   ├── Cronobacter
│       │   │   │   └── Cronobacter sakazakii
│       │   │   ├── Enterobacter
│       │   │   │   ├── Enterobacter asburiae
│       │   │   │   ├── Enterobacter cloacae
│       │   │   │   ├── Enterobacter cloacae complex sp.
│       │   │   │   ├── Enterobacter hormaechei
│       │   │   │   ├── Enterobacter kobei
│       │   │   │   └── Enterobacter roggenkampii
│       │   │   ├── Escherichia
│       │   │   │   ├── Escherichia albertii
│       │   │   │   ├── Escherichia coli
│       │   │   │   ├── Escherichia fergusonii
│       │   │   │   └── Escherichia marmotae
│       │   │   ├── Klebsiella
│       │   │   │   ├── Klebsiella aerogenes
│       │   │   │   ├── Klebsiella grimontii
│       │   │   │   ├── Klebsiella michiganensis
│       │   │   │   ├── Klebsiella oxytoca
│       │   │   │   ├── Klebsiella pneumoniae
│       │   │   │   ├── Klebsiella quasipneumoniae
│       │   │   │   └── Klebsiella variicola
│       │   │   ├── Leclercia
│       │   │   │   └── Leclercia adecarboxylata
│       │   │   ├── Raoultella
│       │   │   │   ├── Raoultella ornithinolytica
│       │   │   │   └── Raoultella planticola
│       │   │   ├── Salmonella
│       │   │   │   ├── Salmonella enterica
│       │   │   │   └── Salmonella sp.
│       │   │   └── Shigella
│       │   │       ├── Shigella boydii
│       │   │       ├── Shigella dysenteriae
│       │   │       ├── Shigella flexneri
│       │   │       └── Shigella sonnei
│       │   ├── Erwiniaceae
│       │   │   ├── Buchnera
│       │   │   │   └── Buchnera aphidicola
│       │   │   ├── Erwinia
│       │   │   │   └── Erwinia amylovora
│       │   │   └── Pantoea
│       │   │       ├── Pantoea agglomerans
│       │   │       └── Pantoea ananatis
│       │   ├── Hafniaceae
│       │   │   └── Edwardsiella
│       │   │       ├── Edwardsiella ictaluri
│       │   │       ├── Edwardsiella piscicida
│       │   │       └── Edwardsiella tarda
│       │   ├── Morganellaceae
│       │   │   ├── Arsenophonus
│       │   │   │   └── Arsenophonus nasoniae
│       │   │   ├── Moellerella
│       │   │   │   └── Moellerella wisconsensis
│       │   │   ├── Morganella
│       │   │   │   └── Morganella morganii
│       │   │   ├── Proteus
│       │   │   │   └── Proteus mirabilis
│       │   │   └── Providencia
│       │   │       ├── Providencia alcalifaciens
│       │   │       ├── Providencia rettgeri
│       │   │       └── Providencia stuartii
│       │   ├── Pectobacteriaceae
│       │   └── Yersiniaceae
│       │       ├── Rahnella
│       │       ├── Serratia
│       │       │   ├── Serratia entomophila
│       │       │   ├── Serratia marcescens
│       │       │   └── Serratia proteamaculans
│       │       └── Yersinia
│       │           ├── Yersinia enterocolitica
│       │           ├── Yersinia pestis
│       │           ├── Yersinia pseudotuberculosis
│       │           └── Yersinia ruckeri
│       ├── Legionellales
│       │   ├── Coxiellaceae
│       │   └── Legionellaceae
│       │       └── Legionella
│       │           ├── Legionella adelaidensis
│       │           └── Legionella pneumophila
│       ├── Methylococcales
│       ├── Moraxellales
│       │   └── Moraxellaceae
│       │       ├── Acinetobacter
│       │       │   ├── Acinetobacter baumannii
│       │       │   ├── Acinetobacter haemolyticus
│       │       │   ├── Acinetobacter indicus
│       │       │   ├── Acinetobacter johnsonii
│       │       │   ├── Acinetobacter lwoffii
│       │       │   ├── Acinetobacter pittii
│       │       │   ├── Acinetobacter schindleri
│       │       │   ├── Acinetobacter seifertii
│       │       │   ├── Acinetobacter ursingii
│       │       │   ├── Acinetobacter variabilis
│       │       │   └── Acinetobacter wuhouensis
│       │       ├── Moraxella
│       │       │   ├── Moraxella bovis
│       │       │   └── Moraxella osloensis
│       │       └── Psychrobacter
│       │           └── Psychrobacter sp.
│       ├── Oceanospirillales
│       ├── Pasteurellales
│       │   └── Pasteurellaceae
│       │       ├── Actinobacillus
│       │       │   └── Actinobacillus pleuropneumoniae
│       │       ├── Glaesserella
│       │       │   └── Glaesserella parasuis
│       │       ├── Haemophilus
│       │       └── Pasteurella
│       │           └── Pasteurella multocida
│       ├── Pseudomonadales
│       │   └── Pseudomonadaceae
│       │       └── Pseudomonas
│       │           ├── Pseudomonas aeruginosa
│       │           ├── Pseudomonas amygdali
│       │           ├── Pseudomonas putida
│       │           ├── Pseudomonas syringae
│       │           └── Pseudomonas syringae group genomosp. 3
│       ├── Thiotrichales
│       │   ├── Francisellaceae
│       │   │   └── Francisella
│       │   └── Piscirickettsiaceae
│       │       └── Piscirickettsia
│       │           └── Piscirickettsia salmonis
│       ├── Vibrionales
│       │   └── Vibrionaceae
│       │       ├── Photobacterium
│       │       │   └── Photobacterium damselae
│       │       └── Vibrio
│       │           ├── Vibrio alginolyticus
│       │           ├── Vibrio campbellii
│       │           ├── Vibrio cholerae
│       │           └── Vibrio parahaemolyticus
│       └── Xanthomonadales
│           └── Xanthomonadaceae
│               ├── Xanthomonas
│               │   ├── Xanthomonas campestris
│               │   ├── Xanthomonas citri
│               │   ├── Xanthomonas hortorum
│               │   ├── Xanthomonas oryzae
│               │   └── Xanthomonas phaseoli
│               └── Xylella
│                   └── Xylella fastidiosa
├── Rhodothermota
│   └── Rhodothermia
│       └── Rhodothermales
│           └── Salinibacteraceae
│               └── Salinibacter
│                   └── Salinibacter ruber
├── Spirochaetota
│   └── Spirochaetia
│       ├── Leptospirales
│       │   └── Leptospiraceae
│       │       └── Leptospira
│       │           ├── Leptospira interrogans
│       │           └── Leptospira noguchii
│       └── Spirochaetales
│           └── Borreliaceae
│               ├── Borrelia
│               │   ├── Borrelia coriaceae
│               │   ├── Borrelia crocidurae
│               │   ├── Borrelia hermsii
│               │   ├── Borrelia miyamotoi
│               │   ├── Borrelia puertoricensis
│               │   └── Borrelia turicatae
│               └── Borreliella
│                   ├── Borreliella bavariensis
│                   ├── Borreliella burgdorferi
│                   ├── Borreliella garinii
│                   └── Borreliella mayonii
└── Thermodesulfobacteriota
    └── Desulfovibrionia
        └── Desulfovibrionales
            └── Desulfovibrionaceae
```

