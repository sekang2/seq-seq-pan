import argparse

from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs=1)

    args = parser.parse_args()

    with open(args.infile[0], "r") as in_file:

        allfiles = in_file.readlines()

        for i in range(0,len(allfiles)):

            fields = allfiles[i].split(";")
            with open(fields[0], "r") as fasta:
                for record in SeqIO.parse(fasta, "fasta"):
                    cur_length = len(record.seq)
                    print("\t".join([str(i+1), record.description, str(cur_length)]))