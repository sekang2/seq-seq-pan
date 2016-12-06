# Pangenome construction


## Prerequisites
This program was implemented in Python and requires Python3.4 or higher.

Software required for running pipeline for set of genomes
* snakemake
* progressiveMauve

---
## Program
### Usage
```
seqseqpan.py  [-h] [-x XMFA_F] -p OUTPUT_P -n OUTPUT_NAME
                [-c CONSENSUS_F] [-m] [-o ORDER]
                [-t {consensus,resolve,realign,xmfa,maf,map,merge,separate}]
                [-i COORD_F] [-l LCB_LENGTH] 

  -h, --help            show this help message and exit
  -x XMFA_F, --xmfa XMFA_F
                        XMFA input file
  -p OUTPUT_P, --output_path OUTPUT_P
                        path to output directory
  -n OUTPUT_NAME, --name OUTPUT_NAME
                        file prefix and sequence header for consensus FASTA /XFMA file
  -c CONSENSUS_F, --consensus CONSENSUS_F
                        consensus FASTA file used in XMFA
  -m, --merge           Merge small blocks to previous or next block in resolve-step.
  -o ORDER, --order ORDER
                        ordering of output (0,1,2,...) [default: 0]
  -t {consensus,resolve,realign,xmfa,maf,map,merge,separate,remove,split}, --task {consensus,resolve,realign,xmfa,maf,map,merge,separate,remove,split}
                        what to do (consensus|resolve|realign|xmfa|map|merge|separate|maf|remove|split)    [default: consensus]
  -i COORD_F, --index COORD_F
                        file with indices to map. First line: source_seq dest_seq[,dest_seq2,...] using "c" or sequence number.
                        Then one coordinate per line. Coordinates are 1-based!
  -l LCB_LENGTH, --length LCB_LENGTH
                        Shorter LCBs will be separated to form genome specific entries.
  -r RM_GENOME, --removegenome RM_GENOME
                        Number of genome to remove (as shown in XMFA header)
```

#### Tasks
Choose task with argument **-t**. Arguments **-p** and **-n** are required for every task.

| task    |description|output|arguments|optional arguments|
|---------|-----------|------|---------|------------------|
|consensus|Create consensus sequence from XMFA file.|2 .FASTA files (with delimiter and without) and 2 .IDX files |-x |-o, -u|
|maf      |Write MAF file from XMFA file.|.MAF file|-x|-o|
|map      |Map positions/coordinates from consensus to sequences, between sequences, ...|.TXT file|-i, -c||
|merge    |Add small LCBs to end or beginning of surrounding LCBs. Stand-alone merging step can only be used with two aligned sequences. |.XMFA file|-x|-o|
|realign  |Realign sequences of LCBs around consecutive gaps, only possible before resolve-step|.XMFA file|-x|-o|
|remove   |Remove a genome from all LCBs in alignment|.XMFA file|-x, -r|-o|
|resolve  |Build alignment of all genomes from .XMFA file with new genome aligned to consensus sequence.|.XMFA file|-x, -c|-o|
|separate |Separate small LCBs to form genome specific entries.|.XMFA file|-x, -l|-o|
|split    |Split LCBs according to chromosome annotation.|.XMFA file|-x|-o|
|xmfa     |Write XMFA file from XMFA file.|.XMFA file|-x|-o|
---

## Pipeline
### Usage
```
snakemake --snakefile run_seqseqpan.Snakemake --config genomefile=genome_list.txt outfilename=TB_example merge=True
```

#### Config

| name        | description |
|-------------|-------------|
| genomefile  |One line per genome with full path and name separated by ';'. |
| outfilename |Prefix for all final output files.|
| merge       |Optional, default = False. Should merging step be included in pipeline?|

---

#### Representation of example
![](representation/example.png)
