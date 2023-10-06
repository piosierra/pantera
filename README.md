![pantera](images/pantera.svg?raw=true "pantera")
# Transposable Elements library generation from a pangenome

A pangenome is a collection of genomes or haplotypes that can be aligned and stored as a variation graph in gfa format. 
:purple_square:**pantera**:purple_square: receives as input a list of gfa files of non overlapping variation graphs and produces a library of transposable elements found to be polymorphic on that pangenome.

## 1- Prepare your gfa files
Use [**pggb**](https://pggb.readthedocs.io/) to create the pangenome from your starting genome sequences. In its most basic form:
1.1- Create a sigle file with the fasta files.
```
cat *.fa combined.fa
```
2.1- Compress and index the file.
```
bgzip -@ 4 combined.fa
samtools faidx combined.fa.gz
```

3.1- Create the pangenome. 
```
pggb \
-i in.fa \       # input file in FASTA format
-o output \      # output directory
-n 9  \          # number of haplotypes
-t 16            # number of threads (defaults to ``getconf _NPROCESSORS_ONLN``)
-p 90 \          # (default) minimum average nucleotide identity for a seed mapping
-s 5000 \        # (default) segment length
```

## 2- Obtain the library from the gfa files
2.1 Create one file with the list of the full paths to the gfa files that will be analyzed.
2.2 Run :purple_square:**pantera**:purple_square:.
```
pantera -g gfas_list -o output_folder
```

## 3- Classify the library obtained
3.1 You can use the TE classifier of your choice. For example, with RepeatClassifier, which is part of the [Dfam TE tools](https://github.com/Dfam-consortium/TETools).
```
RepeatClassifier -consensi pantera_lib.fa
```









