#!/usr/bin/python3

import argparse
import sys
import pdb

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

    if args.task == "map":
        try:
            sparse_align, sparse_consensus = parser.parseConsensusIndex(args.consensus_f)
        except ConsensusFastaIdxFormatError as e:
            print(e.message)
        else:
            try:
                source, dests, coordinates = parser.parseMappingCoordinates(args.coord_f)
            except CoordinatesInputError as e:
                print(e.message)
            else:
                dests.append(source)    
                dests = list(set(dests))
                try:
                    coords_dict = mapper.mapCoordinates(sparse_align, sparse_consensus, source, dests, coordinates)
                except CoordinateOutOfBoundsError as e:
                    print (e.message)
                else:
                    writer.writeMappingCoordinates(source, dests, coords_dict, args.output_p, args.output_name)
    else:
        try:
            
            align = parser.parseXMFA(args.xmfa_f)
                
        except (XMFAHeaderFormatError, LcbInputError) as e:
            print(e.message + "(" + args.xmfa_f + ")")
        else:
            try:
                if args.task == "realign":
                    try: 
                        realign = realigner.realign(align)
                    except ConsensusXMFAInputError as e:
                        print(e.message)
                    else:
                        writer.writeXMFA(realign, args.output_p, args.output_name + "_realign", args.order)
                if args.task == "merge":
                    try:
                        merged = merger.mergeLCBs(align, 1, 2, args.lcb_length)
                        merged = merger.mergeLCBs(merged, 2, 1, args.lcb_length)
                    except ConsensusXMFAInputError as e:
                        print(e.message)
                    else:
                        writer.writeXMFA(merged, args.output_p, args.output_name + "_merge", args.order)
                if args.task == "consensus":
                    
                    writer.writeConsensus(align, args.unambiguous, args.output_p, args.output_name, args.order)
                    
                elif args.task == "separate":
                    if args.lcb_length > 0:
                        separated = separator.separateLCBs(align, args.lcb_length)
                    else:
                        print("Separating LCBs: Nothing do be done -> length is smaller than 1!")
                        separated = align
                        
                    writer.writeXMFA(separated, args.output_p, args.output_name + "_separated", args.order)
                    
                elif args.task == "resolve":
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
                        
                        #pdb.set_trace()
                        resolveblocks_align = resolver.resolveMultiAlignment(align, consensus, consensusGenomeNr=consensusGenomeNr, newGenomeNr=newGenomeNr)
                        #writer.writeXMFA(resolveblocks_align, args.output_p, args.output_name+"_resolvestep", args.order)
                        
                        if args.merge:
                            res_merge = merger.mergeLCBs(resolveblocks_align, consensusGenomeNr=consensusGenomeNr, newGenomeNr=newGenomeNr, blockLength=args.lcb_length)
                            res_merge = merger.mergeLCBs(res_merge, consensusGenomeNr=newGenomeNr, newGenomeNr=consensusGenomeNr, blockLength=args.lcb_length)
                            # realign step necessary in case of consecutive gaps introduced by merging
                            resolveblocks_align = realigner.realign(res_merge)

                            #writer.writeXMFA(resolveblocks_align, args.output_p, args.output_name+"_mergestep", args.order)    
                            
                        reconstruct_align = resolver.reconstructAlignment(resolveblocks_align, consensus, org_align, consensusGenomeNr=consensusGenomeNr, newGenomeNr=newGenomeNr)
                            
                    except (XMFAHeaderFormatError, LcbInputError) as e:
                        print(e.message + "(" + consensus.xmfaFile + ")")
                    except (ConsensusFastaFormatError, ConsensusXMFAInputError, ConsensusGenomeNumberError) as e:
                        print(e.message)
                    else:
                        writer.writeXMFA(reconstruct_align, args.output_p, args.output_name+"_resolve", args.order)
                
                elif args.task == "xmfa":
                    
                    writer.writeXMFA(align, args.output_p, args.output_name, args.order)
                
                elif args.task == "maf":

                    writer.writeMAF(align, args.output_p, args.output_name, args.order)
                    
            except ParameterError as e:
                print('ERROR: Problem with parameter "{0}": Value should be {1}, but was "{2}".'.format(e.parameter, e.rangetext, e.value))
    return(0)
        
if __name__ == '__main__':
               
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file")
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    parser.add_argument("-c", "--consensus", dest="consensus_f", help="consensus FASTA file used in XMFA", required=False)
    parser.add_argument("-u", "--unambiguous", dest="unambiguous", help="Do not use ambigiuous IUPAC code in consensus (random choice instead).", action='store_true')
    parser.add_argument("-m", "--merge", dest="merge", help="Merge small blocks to previous or next block in resolve-step.", action='store_true')
    parser.add_argument("-o", "--order", dest="order", type=int, default=0, help="ordering of output (0,1,2,...) [default: %(default)s]", required=False)
    parser.add_argument("-t", "--task", dest="task", default="consensus", help="what to do (consensus|resolve|realign|xmfa|map|merge|separate|maf) [default: %(default)s]", choices=["consensus", "resolve", "realign", "xmfa", "maf", "map", "merge", "separate"], required=False)
    parser.add_argument("-i", "--index", dest="coord_f", help="file with indices to map. First line: source_seq\tdest_seq[,dest_seq2,...] using \"c\" or sequence number. Then one coordinate per line. Coordinates are 1-based!")
    parser.add_argument("-l", "--length", dest="lcb_length", type=int, help="Shorter LCBs will be separated to form genome specific entries.", required=False, default=10)
    
    args = parser.parse_args()
    
    if args.task == "resolve" and args.consensus_f is None:
         parser.error("Please provide a consensus-sequence file (-c/--consensus) for the \"resolve\"-task (-t/--task).")
    
    if args.task == "map": 
        if args.consensus_f is None:
            parser.error("Please provide a consensus-sequence file (-c/--consensus) for the \"map\"-task (-t/--task).")
        if args.xmfa_f is not None:
            print("WARNING: XMFA file (-x/--xmfa) will not be used for task \"map\" (-t/--task)." , file=sys.stderr)
        if args.coord_f is None:
            parser.error("Please provide a file with indices to map (-i/--index) for task \"map\" (-t/--task).")
    else:
        if args.xmfa_f is None:
            parser.error("Please provide the following arguments: -x/--xmfa")
    
    main()
