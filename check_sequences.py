import argparse
import sys
import pdb


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from seqseqpan.io import Parser


def main():
    global args

    parser = Parser()

    alignment = parser.parse_xmfa(args.xmfa_f)

    genome_dict = {}

    for number, genome in alignment.genomes.items():
        with open(genome.file_path, "r") as fastaf:
            genome_dict[number] = list(SeqIO.parse(fastaf, "fasta"))[0]

    #pdb.set_trace()

    for lcb in alignment.lcbs:
        for entry in lcb.entries:
            seq = str(genome_dict[entry.genome_nr].seq[entry.start-1:entry.end])
            entry_seq = str(entry.sequence.replace("-", ""))
            if entry.strand == "-":
                seq_bio = Seq(entry_seq)
                entry_seq = str(seq_bio.reverse_complement())
            if seq != entry_seq:
                print("Wrong sequence: " + genome_dict[entry.genome_nr].description + "(" + str(entry.genome_nr) + "): " + str(entry.start) + "-" + str(entry.end) + " (" + str(len(entry_seq)) + "," + str(len(seq)) +") "+ str(entry.strand) + "\n" + entry_seq + "\n" + seq + "\n")




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)
    # parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    # parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    args = parser.parse_args()

    main()