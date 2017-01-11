import argparse
import sys
import pdb


from seqseqpan.io import Parser, Writer

def main():
    global args

    parser = Parser()
    writer = Writer()

    alignment = parser.parse_xmfa(args.xmfa_f)
    for lcb in alignment.lcbs:
        lcb.reverse_complement_entries()

    writer.write_xmfa(alignment, args.output_p, args.output_name+"_reversed")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file")
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    args = parser.parse_args()

    main()