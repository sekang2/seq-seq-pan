#!/usr/bin/python3

import argparse
import sys
import pdb

from seqseqpan.io import Parser, Writer
from seqseqpan.modifier import Realigner, Merger, Separator, Remover
from seqseqpan.formatter import Splitter
from seqseqpan.resolver import Resolver
from seqseqpan.exception import *
from seqseqpan.mapper import Mapper


def main():
    global args

    parser = Parser()
    writer = Writer()
    resolver = Resolver()
    realigner = Realigner()
    mapper = Mapper()
    merger = Merger()
    separator = Separator()
    remover = Remover()

    if args.task == "map":
        try:
            sparse_align, sparse_consensus = parser.parse_consensus_index(args.consensus_f)
        except ConsensusFastaIdxFormatError as e:
            print(e.message)
        else:
            try:
                source, destinations_org, coordinates = parser.parse_mapping_coordinates(args.coord_f)
            except CoordinatesInputError as e:
                print(e.message)
            else:
                destinations = list(set(destinations_org)) + [source]
                try:
                    coords_dict = mapper.map_coordinates(sparse_align, sparse_consensus,
                                                         source, destinations, coordinates)
                except CoordinateOutOfBoundsError as e:
                    print(e.message)
                else:
                    writer.write_mapping_coordinates(source, destinations_org, coords_dict, args.output_p, args.output_name)
    else:
        try:
            align = parser.parse_xmfa(args.xmfa_f)

        except (XMFAHeaderFormatError, LcbInputError) as e:
            print(e.message + "(" + args.xmfa_f + ")")
        else:
            try:
                if args.task == "split":

                    chromosome_desc = parser.parse_genome_description(args.genome_desc_f)

                    splitter = Splitter(align, chromosome_desc)
                    split = splitter.split_alignment()

                    writer.write_xmfa(split, args.output_p, args.output_name + "_split", args.order)

                elif args.task == "remove":
                    remove = remover.remove(align, args.rm_genome)
                    merge = remover.merge(remove)

                    writer.write_xmfa(merge, args.output_p, args.output_name + "_removed", args.order)

                elif args.task == "extract":
                    region = args.region

                    region_fields = region.split(":")

                    entries = []
                    for lcb in align.lcbs:
                        entry = lcb.get_entry(int(region_fields[0]))
                        if entry is not None:
                            if entry.strand == "-":
                                entry.reverse_complement()
                            entries.append(entry)

                    sorted_entries = sorted(entries, key=lambda entry: entry.start)

                    sequence = "".join([entry.sequence for entry in sorted_entries])
                    sequence = sequence.replace("-", "")

                    if len(region_fields) > 1:
                        start, end = region_fields[1].split("-")
                        start = int(start)-1
                        sequence = sequence[start:int(end)]


                    writer.write_fasta(region, sequence, args.output_p, args.output_name)

                elif args.task == "realign":
                    try:
                        realign = realigner.realign(align)
                    except ConsensusXMFAInputError as e:
                        print(e.message)
                    else:
                        writer.write_xmfa(realign, args.output_p, args.output_name + "_realign", args.order)
                elif args.task == "merge":
                    try:
                        merged = merger.merge_lcbs(align, 1, 2, args.lcb_length)
                        merged = merger.merge_lcbs(merged, 2, 1, args.lcb_length)
                    except ConsensusXMFAInputError as e:
                        print(e.message)
                    else:
                        writer.write_xmfa(merged, args.output_p, args.output_name + "_merge", args.order)
                elif args.task == "consensus":

                    writer.write_consensus(align, args.output_p, args.output_name, args.order)

                elif args.task == "separate":
                    if args.lcb_length > 0:
                        separated = separator.separate_lcbs(align, args.lcb_length)
                    else:
                        print("Separating LCBs: Nothing do be done -> length is smaller than 1!")
                        separated = align

                    writer.write_xmfa(separated, args.output_p, args.output_name + "_separated", args.order)

                elif args.task == "resolve":
                    try:

                        consensus = parser.parse_block_separated_consensus(args.consensus_f)
                        org_align = parser.parse_xmfa(consensus.xmfa_file)

                        if align.genomes[1].file_path == consensus.fasta_file:
                            consensus_genome_nr = 1
                        elif align.genomes[2].file_path == consensus.fasta_file:
                            consensus_genome_nr = 2
                        else:
                            raise ConsensusGenomeNumberError()

                        new_genome_nr = (1 if consensus_genome_nr == 2 else 2)

                        resolveblocks_align = resolver.resolve_multialignment(align, consensus,
                                                                              consensus_genome_nr=consensus_genome_nr,
                                                                              new_genome_nr=new_genome_nr)

                        if args.merge:
                            res_merge = merger.merge_lcbs(resolveblocks_align, consensus_genome_nr=consensus_genome_nr,
                                                          new_genome_nr=new_genome_nr, block_length=args.lcb_length)
                            res_merge = merger.merge_lcbs(res_merge, consensus_genome_nr=new_genome_nr,
                                                          new_genome_nr=consensus_genome_nr, block_length=args.lcb_length)
                            # realign step necessary in case of consecutive gaps introduced by merging
                            resolveblocks_align = realigner.realign(res_merge)

                        reconstruct_align = resolver.reconstruct_alignment(resolveblocks_align, consensus, org_align,
                                                                           consensus_genome_nr=consensus_genome_nr,
                                                                           new_genome_nr=new_genome_nr)

                    except (XMFAHeaderFormatError, LcbInputError) as e:
                        print(e.message + "(" + consensus.xmfa_file + ")")
                    except (ConsensusFastaFormatError, ConsensusXMFAInputError, ConsensusGenomeNumberError) as e:
                        print(e.message)
                    else:
                        writer.write_xmfa(reconstruct_align, args.output_p, args.output_name + "_resolve", args.order)

                elif args.task == "xmfa":

                    writer.write_xmfa(align, args.output_p, args.output_name, args.order)

                elif args.task == "maf":
                    chromosome_desc = parser.parse_genome_description(args.genome_desc_f)
                    writer.write_maf(align, args.output_p, args.output_name, chromosome_desc)

            except ParameterError as e:
                print('ERROR: Problem with parameter "{0}": Value should be {1}, but was "{2}".'.format(e.parameter,
                                                                                                        e.range_text,
                                                                                                        e.value))
    return 0


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file")
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name",
                        help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    parser.add_argument("-c", "--consensus", dest="consensus_f", help="consensus FASTA file used in XMFA",
                        required=False)
    parser.add_argument("-m", "--merge", dest="merge",
                        help="Merge small blocks to previous or next block in resolve-step.", action='store_true')
    parser.add_argument("-o", "--order", dest="order", type=int, default=0,
                        help="ordering of output (0,1,2,...) [default: %(default)s]", required=False)
    parser.add_argument("-t", "--task", dest="task",
                        help="what to do (consensus|resolve|realign|xmfa|map|merge|separate|maf|remove|split|extract)",
                        choices=["consensus", "resolve", "realign", "xmfa", "maf", "map", "merge", "separate", "remove",
                                 "split", "extract"], required=True)
    parser.add_argument("-i", "--index", dest="coord_f",
                        help="file with indices to map. First line: source_seq\tdest_seq[,dest_seq2,...] using \"c\" or sequence number. Then one coordinate per line. Coordinates are 1-based!")
    parser.add_argument("-l", "--length", dest="lcb_length", type=int,
                        help="Shorter LCBs will be separated to form genome specific entries.", required=False,
                        default=10)
    parser.add_argument("-r", "--removegenome", dest="rm_genome", type=int,
                        help="Number of genome to remove (as shown in XMFA header)", required=False)

    parser.add_argument("-g", "--genome_desc", dest="genome_desc_f", help = "File containing genome description (name/chromosomes) for .MAF file creation and 'split' task.\n"
                                                                    "Format: genome number as in xmfa   name/description    length     (tab-separated)"
                                                                    "Length information is only necessary for FASTA files containing more than one chromosome.\n"
                                                                    "Multiple chromosomes a genome must be listed in the same order as in original FASTA file.\n")
    parser.add_argument("-e", "--extractregion", dest="region", help="Region to extract in the form genome_nr:start-end (one based and inclusive) or only genome_nr for full sequence.")

    args = parser.parse_args()

    if args.task == "resolve" and args.consensus_f is None:
        parser.error("Please provide a consensus-sequence file (-c/--consensus) for the \"resolve\"-task (-t/--task).")

    if args.task == "map":
        if args.consensus_f is None:
            parser.error("Please provide a consensus-sequence file (-c/--consensus) for the \"map\"-task (-t/--task).")
        if args.coord_f is None:
            parser.error("Please provide a file with indices to map (-i/--index) for task \"map\" (-t/--task).")
        if args.xmfa_f is not None:
            print("WARNING: XMFA file (-x/--xmfa) will not be used for task \"map\" (-t/--task).", file=sys.stderr)

    else:
        if args.xmfa_f is None:
            parser.error("Please provide the following arguments: -x/--xmfa")

    if args.task == "remove" and args.rm_genome is None:
        parser.error("Please provide the number of the genome to be removed.")

    if (args.task == "maf" or args.task == "split") and args.genome_desc_f is None:
        parser.error("Please provide the genome description file.")

    if args.task == "extract" and args.region is None:
        parser.error("Please provide the region to extract.")

    main()
