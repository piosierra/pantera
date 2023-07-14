![pantera](images/pantera.svg?raw=true "pantera")
# Transposable Elements library generation from a pangenome

A pangenome is a collection of genomes or haplotypes that can be aligned and stored as a variation graph in gfa format. 
This tool receives as input a list of gfa files of non overlapping variation graphs and produces a library of transposable elements found to be polymorphic on those pangenomes.

DONE:
- Extraction of unique segments
- Consensus sequences for "zones"

TO DO:
- Structure filter
- Final classifier (Optional)
- Final report

## v0.1.4
Fixed defaults.

## v0.1.3
Fixed reverse complement

## v0.1.2
Fast search of segments.
## v0.1.1
- Updated defaults
- Zones definition based only on cluster size
- Fixed number of rounds and behaviour after first
- New consensus function
## v0.0.8
- Added paramether to reduce number of paths to use by size (quantile)
## v0.0.7
- Unified final library
## v0.0.6
- Earlier removal of segments by size
## v0.0.5 
- Improvements on selection of segments
- Correct path added to name
- Zones dynamically created
- More info included in names. It can be used in the postprocessing.
