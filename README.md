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
| task    |description|arguments|optional arguments|
|---------|-----------|---------|------------------|
|consensus|   |   |   |
|maf      |   |   |   |
|map      |   |   |   |
|merge    |   |   |   |
|realign  |   |   |   |
|resolve  |   |   |   |
|separate |   |   |   |
---

## Pipeline
### Usage
snakemake --snakefile run_supergenome.Snakemake --config genomefile=genome_list.txt outfilename=TB_example merge=True

####Config
|name       |description|
|-----------|-----------|
|genomefile |   |
|outfilename|   |
|merge      |   |
