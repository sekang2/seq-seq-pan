#!/usr/bin/python3.4

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
                        merged = merger.mergeLCBs(align, 1, 2)
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
                    
                elif args.task == "split":
                    
                    try:
                        consensus = parser.parseBlockSeparatedConsensus(args.consensus_f)
                        org_align = parser.parseXMFA(consensus.xmfaFile)
                        splitblocks_align = resolver.resolveMultiAlignment(align, consensus, org_align)
                    except (XMFAHeaderFormatError, LcbInputError) as e:
                        print(e.message + "(" + consensus.xmfaFile + ")")
                    except (ConsensusFastaFormatError, ConsensusXMFAInputError, ConsensusGenomeNumberError) as e:
                        print(e.message)
                    else:
                        writer.writeXMFA(splitblocks_align, args.output_p, args.output_name+"_split", args.order)
                
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
    parser.add_argument("-o", "--order", dest="order", type=int, default=0, help="ordering of output (0,1,2,...) [default: %(default)s]", required=False)
    parser.add_argument("-t", "--task", dest="task", default="consensus", help="what to do (consensus|split|realign|xmfa|map|merge|separate|maf) [default: %(default)s]", choices=["consensus", "split", "realign", "xmfa", "maf", "map", "merge", "separate"], required=False)
    parser.add_argument("-i", "--index", dest="coord_f", help="file with indices to map. First line: source_seq\tdest_seq[,dest_seq2,...] using \"c\" or sequence number. Then one coordinate per line.")
    parser.add_argument("-l", "--length", dest="lcb_length", type=int, help="Shorter LCBs will be separated to form genome specific entries.", required=False, default=0)
    
    args = parser.parse_args()
    
    if args.task == "split" and args.consensus_f is None:
         parser.error("Please provide a consensus-sequence file (-c/--consensus) for the \"split\"-task (-t/--task).")
    
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
