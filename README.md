![pantera](images/pantera.png?raw=true "pantera")
# Identification of transposable element families from pangenome polymorphisms

A pangenome is a collection of genomes or haplotypes that can be aligned and stored as a variation graph in gfa format. 
**pantera** receives as input a list of gfa files of non overlapping variation graphs and produces a library of transposable elements found to be polymorphic on that pangenome.

### 0- Installing
Simply download `pantera.R` and make it executable `chmod +x pantera.R` or run with Rscript `Rscript pantera.R` 

### 1- Prepare your gfa files
Use [**pggb**](https://pggb.readthedocs.io/) to create the pangenome from your starting genome sequences. In its most basic form:

1.1- Create a sigle file for each chromosome combining their fasta files.
```
cat *.fa yourspecies.chr1.fa
```
1.2- Compress and index the file.
```
bgzip -@ 4 yourspecies.chr1.fa
samtools faidx yourspecies.chr1.fa.gz
```

1.3- Create a pangenome for each chromosome. `-n` is the number of haplotypes/genomes included. `-t` is the number of threads.
```
pggb -i yourspecies.chr1.fa.gz -o output -n 9 -t 16 
```

### 2- Obtain the library from the gfa files
2.1 Create one file `gfas_list` with the list of the full paths to the gfa files that will be analyzed.

2.2 Run **pantera**, with `-c` as the number of threads.
```
pantera -g gfas_list -c 16 -o output_folder
```

### 2 a- (Optional) Identify structural features in the sequences
```
pantercheck.R pantera_lib.fa
```

### 3- Classify the library obtained
3.1 Use RepeatClassifier, which is part of the [Dfam TE tools](https://github.com/Dfam-consortium/TETools), or your classifier of choice to classify the sequences obtained.
```
RepeatClassifier -consensi pantera_lib.fa
```

### Data example

In the folder `test` there is a gfa example file with the respective outputs that can be used to check if **pantera** works correctly on your system.

***
### Requirements
**pantera** has been tested in Linux with R 4.2.2 to R 4.3.3

**pantera** requires [MAFFT](https://mafft.cbrc.jp/alignment/software/) installed and available in the path.

**pantercheck** requires [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) installed and available in the path.

***
### Citing **pantera**
If you use **pantera** in your work, please cite:

> Sierra, P., Durbin, R. Identification of transposable element families from pangenome polymorphisms.
> *Mobile DNA* **15**, 13 (2024). [doi.org/10.1186/s13100-024-00323-y][doi]












