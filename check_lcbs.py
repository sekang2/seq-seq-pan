import argparse

from seqseqpan.io import Parser, Writer
from seqseqpan.base import *


def main():
    global args

    parser = Parser()
    writer = Writer()

    alignment = parser.parse_xmfa(args.xmfa_f)
    for lcb in alignment.lcbs:
        lcb.entries = sorted(lcb.entries, key=lambda entry: entry.genome_nr)
    new_alignment = Alignment(alignment.xmfa_file)
    for nr, genome in alignment.genomes.items():
        new_alignment.add_genome(genome, nr)
    for order in alignment.genomes:
        if len(alignment.lcbs) > 1:
            alignment.lcbs = alignment.get_sorted_lcbs(order)
            j = 0
            if alignment.lcbs[0].get_entry(order) is not None:
                print("lcb", alignment.lcbs[0].number, "no merge first of genome", order)
                new_alignment.add_lcb(alignment.lcbs[0])
                for lcb in range(1, len(alignment.lcbs)):
                    j += 1
                    if alignment.lcbs[lcb].get_entry(order) is not None:
                        i = 0
                        if len(alignment.lcbs[lcb].entries) == len(alignment.lcbs[lcb-1].entries):
                            strand = alignment.lcbs[lcb].entries[0].strand == alignment.lcbs[lcb-1].entries[0].strand
                            for entry in range(0, len(alignment.lcbs[lcb].entries)):
                                if         alignment.lcbs[lcb].entries[entry].genome_nr != alignment.lcbs[lcb-1].entries[entry].genome_nr\
                                        or alignment.lcbs[lcb].entries[entry].start - alignment.lcbs[lcb-1].entries[entry].end != 1\
                                        or (alignment.lcbs[lcb].entries[entry].strand == alignment.lcbs[lcb-1].entries[entry].strand) != strand:
                                    print("lcb", alignment.lcbs[lcb].number, "no merge1")
                                    new_alignment.add_lcb(alignment.lcbs[lcb])
                                    break
                                else:
                                    i += 1
                            if i == len(alignment.lcbs[lcb].entries):
                                print("lcb", alignment.lcbs[lcb].number, "merge")
                                if not strand:
                                    print("reverse complement lcb", alignment.lcbs[lcb].number)
                                    alignment.lcbs[lcb].reverse_complement_entries()
                                new_alignment.lcbs[-1].length += alignment.lcbs[lcb].length
                                for pos in range(0, len(new_alignment.lcbs[-1].entries)):
                                    new_alignment.lcbs[-1].entries[pos].sequence += alignment.lcbs[lcb].entries[pos].sequence
                                    new_alignment.lcbs[-1].entries[pos].end = alignment.lcbs[lcb].entries[pos].end
                        else:
                            print("lcb", alignment.lcbs[lcb].number, "no merge2")
                            new_alignment.add_lcb(alignment.lcbs[lcb])
                    else:
                        break
            alignment.lcbs[:] = alignment.lcbs[j:len(alignment.lcbs)]
        elif len(alignment.lcbs) == 1 and new_alignment.lcbs[-1].entries[0].genome_nr != alignment.lcbs[0].entries[0].genome_nr:
            print("lcb", alignment.lcbs[0].number, "no merge3")
            new_alignment.add_lcb(alignment.lcbs[0])
            break
        else:
            break

    writer.write_xmfa(new_alignment, args.output_p, args.output_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name",
                        help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    args = parser.parse_args()

    main()
