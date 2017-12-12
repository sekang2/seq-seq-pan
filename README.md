# seq-seq-pan
Build extendable whole genome alignments with a linear represenation optimized for subsequent sequence analysis methods.

## Setup
#### With Conda
To install Conda please refer to the [Conda Installation Guide](https://conda.io/docs/user-guide/install/index.html) and add the [bioconda channel](http://ddocent.com//bioconda/).

Install seq-seq-pan:
```
conda install seq-seq-pan
```

Alternatively, install seq-seq-pan in a separated environment (named ssp here):
```
conda create -n ssp seq-seq-pan
source activate ssp # Activate the environment. To deactivate use "source deactivate"
```
bioconda recipes are also provided as Docker images via [BioContainers](http://biocontainers.pro/).

#### Without Conda
To use seq-seq-pan, the following dependencies have to be satisfied:

* Python3.5 or higher
* Biopython (v1.69) modules: SeqIO, Seq, SeqRecord and pairwise2.
* blat v35
* snakemake (v4.3)
* progressiveMauve (snapshot_2015-02-13)
* Java8

To install seq-seq-pan, run:
```
git clone https://gitlab.com/chrjan/seq-seq-pan.git
cd seq-seq-pan
chmod u+x seq-seq-pan*
```
For easy access, you might want to add the DaisySuite directory to your PATH variable, e.g.
```
export PATH=/path/to/seq-seq-pan/:$PATH
```

## Usage

#### seq-seq-pan-wga
Build pan-genome from set of sequences
```
seq-seq-pan-wga --config genomefile=genome_list.txt outfilename=ssp_example merge=True
```

Add set of sequences to existing pan-genome
```
seq-seq-pan-wga --config genomefile=genome_list_new.txt outfilename=ssp_example_extended merge=True pangenome=ssp_example.xmfa
```

**Configuration**

| name        | description |
|-------------|-------------|
| genomefile  |One line per genome with full path to .FASTA files. |
| outfilename |Prefix for all final output files.|
| merge       |Optional, default = True. Do you want to include the merging steps?|
| pangenome   |Path to exisiting pan-genome (pangenome.XMFA). Accompanying genome description file has to be present in same folder (pangenome_genomedescription.TXT).|
| pmauve      |Optional, path to progressiveMauve executable. If not set, "progressiveMauve" is assumed to be included in $PATH.|
| blat        |Optional, path to blat executable. If not set, "blat" is assumed to be included in $PATH.|

---

#### seq-seq-pan
For full representation of pan-genome and to be able to work with the data structure, only the .XMFA and the genome_description files are needed.
.FASTA files can be discarded as all sequences can be extracted from the pan-genome, the consensus sequence can be constructed from the .XMFA file at any time.

Working with the data structure:

```
seq-seq-pan subcommand ...
```
Get details on subcommand parameters with:
```
seq-seq-pan [subcommand] -h
```

**Subcommands**

Arguments *-p* and *-n* are required for every subcommand.

|subcommand|description|output|arguments|optional arguments|
|----------|-----------|------|---------|------------------|
|blockcountsplit| Split XMFA of 2 genomes into 3 XMFA files: blocks with both genomes and genome-specific blocks for each genome.|3 .XMFA files| -x ||
|consensus      |Create consensus sequence from XMFA file.|2 .FASTA files (with delimiter and without) and 2 .IDX files |-x |-o|
|extract        |Extract sequence for whole genome or genomic interval|.FASTA file|-x, -e||
|join           |Join LCBs from 2 XMFA files, assigning genome_ids as in first XMFA files (-x).|.XMFA file|-x,-y|-o|
|maf            |Write MAF file from XMFA file.|.MAF file|-x, -g||
|map            |Map positions/coordinates from consensus to sequences, between sequences, ...|.TXT file|-i, -c||
|merge          |Add small LCBs to end or beginning of surrounding LCBs. Stand-alone merging step can only be used with two aligned sequences. |.XMFA file|-x|-o, -l|
|realign        |Realign sequences of LCBs around consecutive gaps, only possible before resolve-step.|.XMFA file|-x|-o, --blat|
|reconstruct    |Build alignment of all genomes from .XMFA file with new genome aligned to consensus sequence.|.XMFA file|-x, -c|-o|
|remove         |Remove a genome from all LCBs in alignment.|.XMFA file|-x, -r|-o|
|resolve        |Resolve LCBs stretching over delimiter sequences.|.XMFA file|-x, -c|-o|
|separate       |Separate small LCBs to form genome specific entries.|.XMFA file|-x, -l|-o|
|split          |Split LCBs according to chromosome annotation.|.XMFA file|-x, -g|-o|
|xmfa           |Write XMFA file from XMFA file.|.XMFA file|-x|-o|

---
#### seq-seq-pan-consensus
Build the linear representation of the Pan-genome ("consensus genome").
```
seq-seq-pan-consensus INPUT.xmfa
```

---
#### seq-seq-pan-genomedescription
Generate the genome description file for a set of .FASTA files.
```
seq-seq-pan-genomedescription -i GENOME_LIST -o GENOME_DESC_F
```
