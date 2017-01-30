How to start the example:
-------------------------

unzip seq.zip into seq/.
Check if .fasta paths in "genome_list.txt" are correct.


SIMPLE EXAMPLE:
---------------

Build pan-genome from 6 sequences:
---------------------------------
snakemake --snakefile ../run_seqseqpan.Snakemake --config genomefile=genome_list.txt outfilename=TB_example merge=True



ADVANCED EXAMPLE:
-----------------

Build pan-genome from 4 sequences and add 2 sequences afterwards:
-----------------------------------------------------------------
snakemake --snakefile ../run_seqseqpan.Snakemake --config genomefile=genome_list_4.txt outfilename=TB_example_4 merge=True

snakemake --snakefile ../run_seqseqpan.Snakemake --config genomefile=genome_list_2.txt outfilename=TB_example_6 merge=True pangenome=TB_example_4.xmfa



Add snakemake argument --nt to keep temporary files.

