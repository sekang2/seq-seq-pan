How to start the example:
-------------------------

unzip seq.zip into seq/.
Adjust full path to fasta files in genome_file.txt

Run:
---
snakemake --snakefile /path/to/run_supergenome.Snakemake --config genomefile=/path/to/genome_list.txt outfilename=TB_example merge=True


Add snakemake argument --nt to keep temporary files.