import argparse

from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="genome_list", help= "File with list of /paths/to/files.fasta", required=True)
    parser.add_argument("-o", "--output", dest="genome_desc_f", help="name of output file", required=True)

    args = parser.parse_args()

with open(args.genome_desc_f, "w") as out_file, open(args.genome_list, "r") as in_file:

        allfiles = in_file.readlines()

        for i in range(0,len(allfiles)):

            file = allfiles[i].strip()
            with open(file, "r") as fasta:
                for record in SeqIO.parse(fasta, "fasta"):
                    cur_length = len(record.seq)
                    out_file.write("\t".join([str(i+1), record.description, str(cur_length)]) + "\n")