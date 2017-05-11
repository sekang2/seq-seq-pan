import argparse
import sys
import pdb


from seqseqpan.io import Parser, Writer

def main():
    global args

    parser = Parser()
    writer = Writer()

    alignment = parser.parse_xmfa(args.xmfa_f)

    order = args.order.split(",")
    pdb.set_trace()

    #max_g_len = max([ g.length for g in alignment.genomes.values() ])
    #offset = int("1"+str(max_g_len)) - max_g_len

    for lcb in alignment.lcbs:
        sorted_entries = []
        for genome_nr in order:
            genome_nr = int(genome_nr)
            entry = lcb.get_entry(genome_nr)
            if entry is not None:
                sorted_entries.append(entry)

        lcb.entries = sorted_entries

        if lcb.entries[0].strand == "-":
            lcb.reverse_complement_entries()


    writer.write_xmfa(alignment, args.output_p, args.output_name+"_oriented")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    parser.add_argument("-o", "--order", dest="order", help="comma separated list of order for entries in genomes", required=True)
    args = parser.parse_args()

    main()