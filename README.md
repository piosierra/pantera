![pantera](images/pantera.svg?raw=true "pantera")
# Transposable Elements library generation from a pangenome

A pangenome is a collection of genomes or haplotypes that can be aligned and stored as a variation graph in gfa format. 
**pantera** receives as input a list of gfa files of non overlapping variation graphs and produces a library of transposable elements found to be polymorphic on that pangenome.

### 1- Prepare your gfa files
Use [**pggb**](https://pggb.readthedocs.io/) to create the pangenome from your starting genome sequences. In its most basic form:

1.1- Create a sigle file with the fasta files.
```
cat *.fa combined.fa
```
1.2- Compress and index the file.
```
bgzip -@ 4 combined.fa
samtools faidx combined.fa.gz
```

1.3- Create the pangenome. `-n` is the number of haplotypes/genomes included. `-t` is the number of threads.
```
pggb -i combined.fa.gz -o output -n 9 -t 16 
```

### 2- Obtain the library from the gfa files
2.1 Create one file with the list of the full paths to the gfa files that will be analyzed.
2.2 Run **pantera**, with `-c` as the number of threads.
```
pantera -g gfas_list -c 16 -o output_folder
```

### 3- Classify the library obtained
3.1 You can use the TE classifier of your choice. For example, with RepeatClassifier, which is part of the [Dfam TE tools](https://github.com/Dfam-consortium/TETools).
```
RepeatClassifier -consensi pantera_lib.fa
```









