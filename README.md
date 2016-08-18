# Supergenome construction


## Prerequisites
This program was implemented in Python and requires Python3.4 or higher.

Software required for running pipeline for set of genomes
* snakemake
* progressiveMauve

---
## Program
### Usage
```
supergenome.py  [-h] [-x XMFA_F] -p OUTPUT_P -n OUTPUT_NAME
                [-c CONSENSUS_F] [-u] [-o ORDER]
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
  -u, --unambiguous     Do not use ambigiuous IUPAC code in consensus (random choice instead).
  -o ORDER, --order ORDER
                        ordering of output (0,1,2,...) [default: 0]
  -t {consensus,resolve,realign,xmfa,maf,map,merge,separate}, --task {consensus,resolve,realign,xmfa,maf,map,merge,separate}
                        what to do (consensus|resolve|realign|xmfa|map|merge|separate|maf)    [default: consensus]
  -i COORD_F, --index COORD_F
                        file with indices to map. First line: source_seq dest_seq[,dest_seq2,...] using "c" or sequence number.
                        Then one coordinate per line. Coordinates are 1-based!
  -l LCB_LENGTH, --length LCB_LENGTH
                        Shorter LCBs will be separated to form genome specific entries.
```

#### Tasks
Choose task with argument **-t**. Arguments **-p** and **-n** are required for every task.

| task    |description|output|arguments|optional arguments|
|---------|-----------|------|---------|------------------|
|consensus|Create consensus sequence from XMFA file.|2 .FASTA files (with delimiter and without) and 2 .IDX files |-x |-o, -u|
|maf      |Write MAF file from XMFA file.|.MAF file|-x|-o|
|map      |Map positions/coordinates from consensus to sequences, between sequences, ...|.TXT file|-i, -c||
|merge    |Add small LCBs to end or beginning of surrounding LCBs, only possible before resolve-step.|.XMFA file|-x|-o|
|realign  |Realign sequences of LCBs around consecutive gaps, only possible before resolve-step|.XMFA file|-x|-o|
|resolve  |Build alignment of all genomes from .XMFA file with new genome aligned to consensus sequence.|.XMFA file|-x, -c|-o|
|separate |Separate small LCBs to form genome specific entries.|.XMFA file|-x, -l|-o|
|xmfa     |Write XMFA file from XMFA file.|.XMFA file|-x|-o|
---

## Pipeline
### Usage
```
snakemake --snakefile run_supergenome.Snakemake --config genomefile=genome_list.txt outfilename=TB_example merge=True
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