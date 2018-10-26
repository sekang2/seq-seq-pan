#!/usr/bin/env python3

import argparse

from seqseqpan.io import Parser, Writer, Processor
from seqseqpan.modifier import Realigner, Merger, Separator, Remover, SingletonAligner
from seqseqpan.formatter import Splitter
from seqseqpan.resolver import Resolver
from seqseqpan.exception import *
from seqseqpan.mapper import Mapper

def parse_xmfa(xmfa_f):
    try:
        align = parser.parse_xmfa(xmfa_f)
    except (XMFAHeaderFormatError, LcbInputError) as e:
        print(e.message + "(" + xmfa_f + ")")
    return align


def parse_consensus(align, consensus_f):
    try:
        consensus = parser.parse_block_separated_consensus(consensus_f)

        if align.genomes[1].file_path == consensus.fasta_file:
            consensus_genome_nr = 1
        elif align.genomes[2].file_path == consensus.fasta_file:
            consensus_genome_nr = 2
        else:
            raise ConsensusGenomeNumberError()

        new_genome_nr = (1 if consensus_genome_nr == 2 else 2)

    except (ConsensusFastaFormatError, ConsensusXMFAInputError, ConsensusGenomeNumberError) as e:
        print(e.message)

    return consensus, consensus_genome_nr, new_genome_nr


# functions for tasks

def blockcountsplit_task(args):
    singleton_aligner = SingletonAligner()
    align = parse_xmfa(args.xmfa_f)
    pairblocks_alignment, consensus_singleton_alignment, new_singleton_alignment = \
        singleton_aligner.genome_count_split(align)

    writer.write_xmfa(consensus_singleton_alignment, args.output_p,
                      args.output_name + "_1_single",
                      order=0, check_invalid=False)
    writer.write_xmfa(new_singleton_alignment, args.output_p,
                      args.output_name + "_2_single",
                      order=0, check_invalid=False)
    writer.write_xmfa(pairblocks_alignment, args.output_p,
                      args.output_name + "_pairblocks",
                      order=0, check_invalid=False)


def consensus_task(args):
    align = parse_xmfa(args.xmfa_f)
    writer.write_consensus(align, args.output_p, args.output_name, args.order)


def extract_task(args):
    align = parse_xmfa(args.xmfa_f)

    region = args.region

    chromosome_desc = parser.parse_genome_description(args.genome_desc_f)
    splitter = Splitter(align, chromosome_desc)

    chromosomes, chunks = splitter.split_sequence(region)

    writer.write_fasta(chromosomes, chunks, args.output_p, args.output_name)


def join_task(args):
    singleton_aligner = SingletonAligner()
    align = parse_xmfa(args.xmfa_f)
    try:
        align_2 = parser.parse_xmfa(args.xmfa_f_2)
    except (XMFAHeaderFormatError, LcbInputError) as e:
        print(e.message + "(" + args.xmfa_f_2 + ")")
    else:
        joined_alignment = singleton_aligner.join(align, align_2)
        writer.write_xmfa(joined_alignment, args.output_p, args.output_name + "_joined", args.order,
                          check_invalid=(not args.quiet))


def maf_task(args):
    align = parse_xmfa(args.xmfa_f)
    chromosome_desc = parser.parse_genome_description(args.genome_desc_f)
    writer.write_maf(align, args.output_p, args.output_name, chromosome_desc, check_invalid=True)


def map_task(args):
    mapper = Mapper()
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


def merge_task(args):
    merger = Merger()
    align = parse_xmfa(args.xmfa_f)
    try:
        merged = merger.merge_lcbs(align, 1, 2, args.lcb_length)
        merged = merger.merge_lcbs(merged, 2, 1, args.lcb_length)
    except ConsensusXMFAInputError as e:
        print(e.message)
    else:
        writer.write_xmfa(merged, args.output_p, args.output_name + "_merge", args.order, check_invalid=(not args.quiet))


def realign_task(args):
    align = parse_xmfa(args.xmfa_f)
    realigner = Realigner()
    processor = Processor(args.output_p, blat=args.blat_path)
    try:
        realign = realigner.realign(align, processor)
    except ConsensusXMFAInputError as e:
        print(e.message)
    else:
        writer.write_xmfa(realign, args.output_p, args.output_name + "_realign", args.order, check_invalid=(not args.quiet))


def reconstruct_task(args):
    resolver = Resolver()
    align = parse_xmfa(args.xmfa_f)
    consensus, consensus_genome_nr, new_genome_nr = parse_consensus(align, args.consensus_f)
    try:
        org_align = parser.parse_xmfa(consensus.xmfa_file)

        reconstruct_align = resolver.reconstruct_alignment(align, consensus, org_align,
                                                           consensus_genome_nr=consensus_genome_nr,
                                                           new_genome_nr=new_genome_nr)
    except (XMFAHeaderFormatError, LcbInputError) as e:
        print(e.message + "(" + consensus.xmfa_file + ")")
    else:
        writer.write_xmfa(reconstruct_align, args.output_p, args.output_name + "_reconstruct",
                          args.order, check_invalid=(not args.quiet))


def remove_task(args):
    remover = Remover()
    align = parse_xmfa(args.xmfa_f)

    remove = remover.remove(align, args.rm_genome)
    merged = remover.merge(remove)

    writer.write_xmfa(merged, args.output_p, args.output_name + "_removed", args.order, check_invalid=(not args.quiet))


def resolve_task(args):
    resolver = Resolver()
    align = parse_xmfa(args.xmfa_f)

    consensus, consensus_genome_nr, new_genome_nr = parse_consensus(align, args.consensus_f)

    resolveblocks_align = resolver.resolve_multialignment(align, consensus,
                                                          consensus_genome_nr=consensus_genome_nr,
                                                          new_genome_nr=new_genome_nr)
    writer.write_xmfa(resolveblocks_align, args.output_p, args.output_name + "_resolve",
                      args.order, check_invalid=False)


def separate_task(args):
    separator = Separator()
    align = parse_xmfa(args.xmfa_f)
    if args.lcb_length > 0:
        separated = separator.separate_lcbs(align, args.lcb_length)
    else:
        print("Separating LCBs: Nothing do be done -> length is smaller than 1!")
        separated = align

    writer.write_xmfa(separated, args.output_p, args.output_name + "_separated", args.order, check_invalid=(not args.quiet))


def split_task(args):
    align = parse_xmfa(args.xmfa_f)
    chromosome_desc = parser.parse_genome_description(args.genome_desc_f)

    splitter = Splitter(align, chromosome_desc)
    split = splitter.split_alignment()

    writer.write_xmfa(split, args.output_p, args.output_name + "_split", args.order, check_invalid=(not args.quiet))


def xmfa_task(args):
    align = parse_xmfa(args.xmfa_f)
    writer.write_xmfa(align, args.output_p, args.output_name, args.order, check_invalid=(not args.quiet))


if __name__ == '__main__':
    # top parser with always needed arguments
    args_parser = argparse.ArgumentParser(description="seq-seq-pan - Working with the Pan-genome")

    required = args_parser.add_argument_group('required arguments')

    # these arguments are required for all parsers - but not shown in help for subparser
    # they are therefore moved to a "parent parser" and included in 'parent' for all subparsers
    top_parser = argparse.ArgumentParser(add_help=False)
    top_parser.add_argument("--quiet", dest="quiet", help="Suppress warnings.", action='store_true')
    top_required = top_parser.add_argument_group('required arguments')

    top_required.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    top_required.add_argument("-n", "--name", dest="output_name",
                        help="File prefix and sequence header for output FASTA / XFMA file", required=True)


    # parent parser for common arguments
    xmfa_file_parser = argparse.ArgumentParser(add_help=False)
    xmfa_required = xmfa_file_parser.add_argument_group("required arguments")
    xmfa_required.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)

    order_parser = argparse.ArgumentParser(add_help=False)
    order_parser.add_argument("-o", "--order", dest="order", type=int, default=0,
                        help="Ordering of blocks in XMFA/FASTA output (0,1,2,...) [default: %(default)s]", required=False)

    genomedesc_parser = argparse.ArgumentParser(add_help=False)
    genomedesc_required = genomedesc_parser.add_argument_group("required arguments")
    genomedesc_required.add_argument("-g", "--genome_desc", dest="genome_desc_f",
                        help="File containing genome description (name/chromosomes) for .MAF file creation and 'split' task. \n"
                             "Format: 'genome id<TAB>name<TAB>length' \n"
                             "Length information is only necessary for FASTA files containing more than one chromosome. \n"
                             "Multiple chromosomes of a genome must be listed in the same order as in original FASTA file.\n",
                        required=True)

    consensus_file_parser = argparse.ArgumentParser(add_help=False)
    consensus_file_required = consensus_file_parser.add_argument_group("required arguments")
    consensus_file_required.add_argument("-c", "--consensus", dest="consensus_f", help="consensus FASTA file used in XMFA",
                                       required=True)

    block_length_parser = argparse.ArgumentParser(add_help=False)
    block_length_parser.add_argument("-l", "--length", dest="lcb_length", type=int,
                                 help="Length of \"small LCB\". [default: %(default)s]", required=False, default=10)


    # subparser for each task
    subparsers = args_parser.add_subparsers(title="subcommands", dest="subcommand", metavar="subcommand",
                                       description="")
    subparsers.required = True

    # blockcountsplit
    subcommand_text = "Split XMFA of 2 genomes into 3 XMFA files: blocks with both genomes and genome-specific blocks for each genome."
    blockcountsplit_parser = subparsers.add_parser("blockcountsplit", help=subcommand_text, description=subcommand_text,
                                                   parents=[top_parser, xmfa_file_parser])
    blockcountsplit_parser.set_defaults(func=blockcountsplit_task)

    # consensus
    #subcommand_text = "Create consensus sequence from XMFA file."
    #consensus_parser = subparsers.add_parser("consensus", help=subcommand_text, description=subcommand_text,
    #                                         parents=[top_parser, xmfa_file_parser, order_parser])
    #consensus_parser.set_defaults(func=consensus_task)

    # extract
    subcommand_text = "Extract sequence for whole genome or genomic interval."
    extract_parser = subparsers.add_parser("extract", help=subcommand_text, description=subcommand_text,
                                           parents=[top_parser, xmfa_file_parser, genomedesc_parser])
    extract_parser._action_groups[2].add_argument("-e", "--extractregion", dest="region",
                        help="Region to extract in the form genome_nr:start-end (one based and inclusive) "
                             "or only genome_nr for full sequence.", required=True)
    extract_parser.set_defaults(func=extract_task)

    # join
    subcommand_text = "Join LCBs from 2 XMFA files, assigning genome_ids as in first XMFA file (-x)."
    join_parser = subparsers.add_parser("join", help=subcommand_text, description=subcommand_text,
                                        parents=[top_parser, xmfa_file_parser, order_parser])
    join_parser._action_groups[2].add_argument("-y", "--xmfa_two", dest="xmfa_f_2",
                                               help="XMFA file to be joined with input file.", required=True)
    join_parser.set_defaults(func=join_task)

    # maf
    subcommand_text = "Write MAF file from XMFA file."
    maf_parser = subparsers.add_parser("maf", help=subcommand_text, description=subcommand_text,
                                       parents=[top_parser, xmfa_file_parser, genomedesc_parser])
    maf_parser.set_defaults(func=maf_task)

    # map
    subcommand_text = "Map positions/coordinates from consensus to sequences, between sequences, ..."
    map_parser = subparsers.add_parser("map", help=subcommand_text, description=subcommand_text,
                                       parents=[top_parser, consensus_file_parser])
    map_parser._action_groups[2].add_argument("-i", "--index", dest="coord_f",
                        help="File with indices to map. First line: source_seq<TAB>dest_seq[,dest_seq2,...] using \"c\" or "
                             "sequence number. Then one coordinate per line. Coordinates are 1-based!", required=True)
    map_parser.set_defaults(func=map_task)

    # merge
    subcommand_text = "Add small LCBs to end or beginning of surrounding LCBs. Can only be used with two aligned sequences."
    merge_parser = subparsers.add_parser("merge", help=subcommand_text, description=subcommand_text,
                                         parents=[top_parser, xmfa_file_parser, order_parser, block_length_parser])
    merge_parser.set_defaults(func=merge_task)

    # realign
    subcommand_text = "Realign sequences of LCBs around consecutive gaps. Can only be used with two aligned sequences."
    realign_parser = subparsers.add_parser("realign", help=subcommand_text, description=subcommand_text,
                                           parents=[top_parser, xmfa_file_parser, order_parser])
    realign_parser.add_argument("--blat", dest="blat_path", help="Path to blat binary if not in $PATH.", default="blat")
    realign_parser.set_defaults(func=realign_task)

    # reconstruct
    subcommand_text = "Build alignment of all genomes from .XMFA file with new genome aligned to consensus sequence."
    reconstruct_parser = subparsers.add_parser("reconstruct", help=subcommand_text, description=subcommand_text,
                                               parents=[top_parser, xmfa_file_parser, consensus_file_parser, order_parser])
    reconstruct_parser.set_defaults(func=reconstruct_task)

    # remove
    subcommand_text = "Remove a genome from all LCBs in alignment."
    remove_parser = subparsers.add_parser("remove", help=subcommand_text, description=subcommand_text,
                                          parents=[top_parser, xmfa_file_parser, order_parser])
    remove_parser._action_groups[2].add_argument("-r", "--removegenome", dest="rm_genome", type=int,
                        help="Number of genome to remove (as shown in XMFA header)", required=True)
    remove_parser.set_defaults(func=remove_task)

    # resolve
    subcommand_text = "Resolve LCBs stretching over delimiter sequences."
    resolve_parser = subparsers.add_parser("resolve", help=subcommand_text, description=subcommand_text,
                                           parents=[top_parser, xmfa_file_parser, consensus_file_parser, order_parser])
    resolve_parser.set_defaults(func=resolve_task)

    # separate
    subcommand_text = "Separate small multi-sequence LCBs to form genome specific entries."
    separate_parser = subparsers.add_parser("separate", help=subcommand_text, description=subcommand_text,
                                            parents=[top_parser, xmfa_file_parser, order_parser, block_length_parser])
    separate_parser.set_defaults(func=separate_task)

    # split
    subcommand_text = "Split LCBs according to chromosome annotation."
    split_parser = subparsers.add_parser("split", help=subcommand_text, description=subcommand_text,
                                         parents=[top_parser, xmfa_file_parser, genomedesc_parser, order_parser])
    split_parser.set_defaults(func=split_task)

    # xmfa
    subcommand_text = "Write XMFA file from XMFA file (for reordering or checking validity)."
    xmfa_parser = subparsers.add_parser("xmfa", help=subcommand_text, description=subcommand_text,
                                        parents=[top_parser, xmfa_file_parser, order_parser])
    xmfa_parser.set_defaults(func=xmfa_task)

    # global seqseqpan parser and writer
    parser = Parser()
    writer = Writer()

    # parse all arguments and call subparser function
    all_args = args_parser.parse_args()
    all_args.func(all_args)

