# HOTSPOT
HOTSPOT is a hierarchical host prediction tool for plasmid contigs using Transformer.

HOTSPOT can learn the gene content features of the input DNA sequences based on deep learning classifiers. Using these features, HOTSPOT will predict the taxonomic lineage information of the input. There are 115 classifier models in all, each corresponding to a taxonomy tree node. To use HOTSPOT, you only need to input contigs of plasmid origin into the program.

# Required Dependencies




# Prepare the database
The database of HOTSPOT will occupy a large storage space (45G) because of the magnitude of the full plasmid taxonomy tree. To 


# Format of the output file
1. The format of the input file should be a fasta file that contains one or more plasmid contigs (complete plasmids are also ok). If an input contig is too short (e.g., <1000 bp), it's possible that **Prodigal** can predict no protein. In this case, HOTSPOT will not output any result for this contig.
2. The output file format is a csv file containing the prediction taxonomic lineage information of input plasmid contigs. Each row represents the prediction result for one input contig. **Contig** is the contig accession number from the input.
| Contig | phylum | class | order | family | genus | species |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| CP090713.1  | Proteobacteria  | Betaproteobacteria  | Burkholderiales  | Burkholderiaceae  | Burkholderia  | Burkholderia multivorans  |
| NZ_CP050042.1  | Proteobacteria  | Gammaproteobacteria  | Enterobacterales  | Enterobacteriaceae  | Escherichia  | Escherichia coli  | 
| NZ_CP083619.1  | Firmicutes  | Clostridia  | Eubacteriales  | Peptostreptococcaceae  | Clostridioides  | Clostridioides difficile  |
| NZ_CP083659.1  | Proteobacteria  | Gammaproteobacteria  | Moraxellales  | Moraxellaceae  | Acinetobacter  | Acinetobacter variabilis  |
| Z22927.1  | Actinobacteria  | Actinomycetia  | Corynebacteriales  | Corynebacteriaceae  | Corynebacterium  | Corynebacterium glutamicum  |



# Contact
If you have any questions, please email us:
  
>yongxinji2-c@my.cityu.edu.hk (Ji Yongxin)
  
>jyshang2-c@my.cityu.edu.hk (Shang Jiayu)
  
>xubotang2-c@my.cityu.edu.hk (Tang Xubo)
