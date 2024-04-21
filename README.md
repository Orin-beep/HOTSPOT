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
* For GPU version PyTorch: search [PyTorch](https://pytorch.org/get-started/previous-versions/) to find the correct CUDA version according to your computer


## Prepare the environment
After cloning this repository (```git clone https://github.com/Orin-beep/HOTSPOT```), you can use Anaconda to install ```environment.yaml```. This will install all packages you need in GPU mode (make sure you have installed CUDA on your system to use the GPU version; otherwise, HOTSPOT will run in CPU mode). The installing command is: 
```
cd HOTSPOT/
conda env create -f environment.yaml -n hotspot
conda activate hotspot
```
If Anaconda fails to work, you can prepare the environment by individually installing the packages listed in the __Dependencies__ section.

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


## Format of the output file
The results will be saved in a TSV file containing the predicted host lineages from phylum to species level. Each row corresponds to an input plasmid sequence. Examples:

| Contig | Phylum | Class | Order | Family | Genus | Species |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| NZ_CP050042.1  | Pseudomonadota  | Gammaproteobacteria  | Enterobacterales  | Enterobacteriaceae  | Escherichia  | -  | 
| NZ_CP083619.1  | Bacillota  | Clostridia  | Eubacteriales  | Peptostreptococcaceae  | Clostridioides  | Clostridioides difficile  |
| NZ_CP083659.1  | Pseudomonadota  | Gammaproteobacteria  | Moraxellales  | Moraxellaceae  | Acinetobacter  | Acinetobacter variabilis  |
| Z22927.1  | Actinomycetota  | Actinomycetes  | Corynebacteriales  | Corynebacteriaceae  | Corynebacterium  | Corynebacterium glutamicum  |

The dash '-' indicates that HOTSPOT cannot provide accurate predictions for the corresponding input at this taxonomic level.


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

