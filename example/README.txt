How to start the example:
-------------------------

Create directory seq/ in current directory and unzip seq.zip into seq/.


SIMPLE EXAMPLE:
---------------

Build pan-genome from 6 sequences:
---------------------------------
../seq-seq-pan-wga --config genomefile=genome_list.txt outfilename=ssp_example



ADVANCED EXAMPLE:
-----------------

Build pan-genome from 4 sequences and add 2 sequences afterwards:
-----------------------------------------------------------------
../seq-seq-pan-wga --config genomefile=genome_list_4.txt outfilename=ssp_example_4 merge=True

../seq-seq-pan-wga --config genomefile=genome_list_2.txt outfilename=ssp_example_6 merge=True pangenome=ssp_example_4.xmfa



Add snakemake argument --nt to keep temporary files.

