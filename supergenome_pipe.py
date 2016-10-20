#!/usr/bin/python3

import argparse
import sys

from supergenome.io import Parser, Writer
from supergenome.modifier import Realigner, Merger, Separator
from supergenome.resolver import Resolver
from supergenome.exception import *
from supergenome.mapper import Mapper

def main():

    global args
    
    parser = Parser()
    writer = Writer()
    resolver = Resolver()
    realigner = Realigner()
    mapper = Mapper()
    merger = Merger()
    separator = Separator() 
    
    try:
            
        align = parser.parseXMFA(args.xmfa_f)
            
    except (XMFAHeaderFormatError, LcbInputError) as e:
        print(e.message + "(" + args.xmfa_f + ")")
    else:
    
        final_align = align
        
        try:
            if args.merge:
                final_align = merger.mergeLCBs(final_align, 1, 2, args.lcb_length)
                final_align = merger.mergeLCBs(final_align, 2, 1, args.lcb_length)
                
            final_align = realigner.realign(final_align)
        except ConsensusXMFAInputError as e:
            print(e.message)
        else:
        
            if args.consensus_f is not None:
        
                try:
                    consensus = parser.parseBlockSeparatedConsensus(args.consensus_f)
                    org_align = parser.parseXMFA(consensus.xmfaFile)
                    
                    if align.genomes[1].filepath == consensus.fastaFile:
                        consensusGenomeNr = 1
                    elif align.genomes[2].filepath == consensus.fastaFile:
                        consensusGenomeNr = 2
                    else:
                        raise ConsensusGenomeNumberError()
                    
                    newGenomeNr = (1 if consensusGenomeNr == 2 else 2)
                    
                    resolveblocks_align = resolver.resolveMultiAlignment(final_align, consensus, consensusGenomeNr=consensusGenomeNr, newGenomeNr=newGenomeNr)
                    
                    if args.merge:
                        res_merge = merger.mergeLCBs(resolveblocks_align, consensusGenomeNr=consensusGenomeNr, newGenomeNr=newGenomeNr, blockLength=args.lcb_length)
                        # realign step necessary in case of consecutive gaps introduced by merging
                        resolveblocks_align = realigner.realign(res_merge)

                    final_align = resolver.reconstructAlignment(resolveblocks_align, consensus, org_align, consensusGenomeNr=consensusGenomeNr, newGenomeNr=newGenomeNr)
                        
                except (XMFAHeaderFormatError, LcbInputError) as e:
                    print(e.message + "(" + consensus.xmfaFile + ")")
                except (ConsensusFastaFormatError, ConsensusXMFAInputError, ConsensusGenomeNumberError) as e:
                    print(e.message)
            
            writer.writeXMFA(final_align, args.output_p, args.output_name+"_ready")
    return 0 
    
    
if __name__ == '__main__':
               
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    parser.add_argument("-c", "--consensus", dest="consensus_f", help="consensus FASTA file used in XMFA", required=False)
    parser.add_argument("-m", "--merge", dest="merge", help="Merge small blocks to previous or next block in resolve-step.", action='store_true')
    parser.add_argument("-l", "--length", dest="lcb_length", type=int, help="Shorter LCBs will be separated to form genome specific entries.", required=False, default=10)

    args = parser.parse_args()
    
    main()