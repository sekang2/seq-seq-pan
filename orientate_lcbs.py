import argparse
import sys
import pdb


from seqseqpan.io import Parser, Writer

def main():
    global args

    parser = Parser()
    writer = Writer()

    alignment = parser.parse_xmfa(args.xmfa_f)

    max_g_len = max([ g.length for g in alignment.genomes.values() ])
    offset = int("1"+str(max_g_len)) - max_g_len

    for lcb in alignment.lcbs:
        lcb.entries = sorted(lcb.entries, key=lambda entry, offset=offset:
                                     entry.start + offset*entry.genome_nr)
        if lcb.entries[0].strand == "-":
            lcb.reverse_complement_entries()


    writer.write_xmfa(alignment, args.output_p, args.output_name+"_oriented")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    args = parser.parse_args()

    main()