import argparse

from Bio import SeqIO


def extract_description(genome_list, genome_desc_f, add_f=None):
    with open(genome_desc_f, "w") as out_file, open(genome_list, "r") as in_file:

        offset = 0

        if add_f is not None:
            with open(add_f) as add:
                all_lines = add.readlines()
                offset = len(all_lines)
                out_file.writelines(all_lines)

        allfiles = in_file.readlines()

        for i in range(0, len(allfiles)):

            file = allfiles[i].strip()
            with open(file, "r") as fasta:
                for record in SeqIO.parse(fasta, "fasta"):
                    cur_length = len(record.seq)
                    out_file.write("\t".join([str(i + 1 + offset), record.description, str(cur_length)]) + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="genome_list", help= "File with list of /paths/to/files.fasta", required=True)
    parser.add_argument("-o", "--output", dest="genome_desc_f", help="name of output file", required=True)
    parser.add_argument("-a", "--add", dest="add_f",  help="Add new genome description to this file.")

    args = parser.parse_args()

    extract_description(args.genome_list, args.genome_desc_f, args.add_f)
