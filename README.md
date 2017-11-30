# Pan-genome construction


## Prerequisites
This program was implemented in Python and requires Python3.4 or higher.<br/>

It depends on the following Biopython modules: SeqIO, Seq, SeqRecord and pairwise2.<br/>
It depends on the following software: blat <br\>

<br/>  
Software required for running pipeline for set of genomes
* snakemake
* progressiveMauve
* Java8


## Pipeline
### Usage
#### Build pan-genome from set of sequences
```
snakemake --snakefile run_seqseqpan.Snakemake --config genomefile=genome_list.txt outfilename=TB_example merge=True
```

#### Add set of sequences to existing pan-genome
```
snakemake --snakefile run_seqseqpan.Snakemake --config genomefile=genome_list_new.txt outfilename=TB_example_extended merge=True pangenome=TB_example.xmfa
```

#### Config

| name        | description |
|-------------|-------------|
| genomefile  |One line per genome with full path to .FASTA files. |
| outfilename |Prefix for all final output files.|
| merge       |Optional, default = True. Do you want to include the merging steps?|
| pangenome   |Path to exisiting pan-genome (pangenome.XMFA). Accompanying genome description file has to be present in same folder (pangenome_genomedescription.TXT).|
| pmauve      |Optional, path to progressiveMauve executable. If not set, "progressiveMauve" is assumed to be included in $PATH.|
| blat        |Optional, path to blat executable. If not set, "blat" is assumed to be included in $PATH.|

---


## Pan-genome data structure
For full representation of pan-genome and to be able to work with the data structure, only the .XMFA and the genome_description files are needed.

.FASTA files can be discarded as all sequences can be extracted from the pan-genome, the consensus sequence can be constructed from the .XMFA file at any time.


### Working with the data structure

```
usage: seqseqpan.py [-h] subcommand ...
```

#### Subcommands
Arguments **-p** and **-n** are required for every subcommand.

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

Get details on subcommand parameters with:
```
seqseqpan.py [subcommand] -h
```


#### Additional scripts
##### Genome Description File
Use the genome_description.py script to generate the genome description file for a set of .FASTA files.
```
python3.4 genomedescription.py -i GENOME_LIST -o GENOME_DESC_F
```
