#!/usr/bin/python3

import argparse
import sys

from seqseqpan.io import Parser, Writer
from seqseqpan.modifier import Realigner, Merger
from seqseqpan.resolver import Resolver
from seqseqpan.exception import *


def main():
    global args

    parser = Parser()
    writer = Writer()
    resolver = Resolver()
    realigner = Realigner()
    merger = Merger()

    try:

        align = parser.parse_xmfa(args.xmfa_f)

    except (XMFAHeaderFormatError, LcbInputError) as e:
        print(e.message + "(" + args.xmfa_f + ")")
    else:

        final_align = align

        try:
            if args.merge:
                final_align = merger.merge_lcbs(final_align, 1, 2, args.lcb_length)
                final_align = merger.merge_lcbs(final_align, 2, 1, args.lcb_length)

            final_align = realigner.realign(final_align)
        except ConsensusXMFAInputError as e:
            print(e.message)
        else:

            if args.consensus_f is not None:

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

                    resolve_blocks_align = resolver.resolve_multialignment(final_align, consensus,
                                                                           consensus_genome_nr=consensus_genome_nr,
                                                                           new_genome_nr=new_genome_nr)

                    if args.merge:
                        res_merge = merger.merge_lcbs(resolve_blocks_align, consensus_genome_nr=consensus_genome_nr,
                                                      new_genome_nr=new_genome_nr, block_length=args.lcb_length)
                        # realign step necessary in case of consecutive gaps introduced by merging
                        resolve_blocks_align = realigner.realign(res_merge)

                    final_align = resolver.reconstruct_alignment(resolve_blocks_align, consensus, org_align,
                                                                 consensus_genome_nr=consensus_genome_nr,
                                                                 new_genome_nr=new_genome_nr)

                except (XMFAHeaderFormatError, LcbInputError) as e:
                    print(e.message + "(" + consensus.xmfa_file + ")")
                except (ConsensusFastaFormatError, ConsensusXMFAInputError, ConsensusGenomeNumberError) as e:
                    print(e.message)

            writer.write_xmfa(final_align, args.output_p, args.output_name + "_ready")
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name",
                        help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    parser.add_argument("-c", "--consensus", dest="consensus_f",
                        help="consensus FASTA file used in XMFA", required=False)
    parser.add_argument("-m", "--merge", dest="merge",
                        help="Merge small blocks to previous or next block in resolve-step.", action='store_true')
    parser.add_argument("-l", "--length", dest="lcb_length", type=int,
                        help="Shorter LCBs will be separated to form genome specific entries.",
                        required=False, default=10)

    args = parser.parse_args()

    main()
