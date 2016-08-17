import sys
import os
import re
import random

from Bio import SeqIO
from Bio.Seq import Seq

from supergenome.exception import *
from supergenome.constants import BLOCK_DELIMITER, RANDOM_SEED

class Genome:

    def __init__(self, filepath, format, entry=-1):
        self.filepath = os.path.abspath(filepath)
        self.format = format
        self.entry = int(entry)
        self.chromosomes = None
        
    def readChromosomes(self):
        self.chromosomes = {}
        start = 1
        with open(self.filepath, "r") as fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                cur_length = len(record.seq)
                self.chromosomes[int(start)] = {"desc": record.description, "length": cur_length}
                start += cur_length
                
                
class Alignment:

    def __init__(self, xmfaFile):
        self.xmfaFile = os.path.abspath(xmfaFile)
        self.LCBs = []
        self.genomes = {}
        
        
    def addGenome(self, genome, number):
        self.genomes[int(number)] = genome
        
        
    def addLCBEntries(self, entries):
        number = len(self.LCBs) + 1
        lcb = LCB(number)
        lcb.addEntries(entries)
        self.LCBs.append(lcb)
        
        
    def addLCB(self, lcb):
        number = len(self.LCBs) + 1
        lcb.number = number
        self.LCBs.append(lcb)
        
        
    def getConsensus(self, unambiguous, order=0):
        sortedLCBs = self.getSortedLCBs(order)    
        delim = BLOCK_DELIMITER
        
        consensusseqs = [lcb.consensusSequence(unambiguous) for lcb in sortedLCBs]
        
        #store beginning of each block calculated using lengths of joined LCBs
        lcblengths = [lcb.length for lcb in sortedLCBs]
        lcblengths = [0] + lcblengths
        for i in range(1,len(lcblengths)):
            lcblengths[i] = lcblengths[i-1] + lcblengths[i] + len(delim)
        
        consseq = delim.join( consensusseqs )
        
        return consseq, lcblengths[:-1]
        

    def getSortedLCBs(self, order):
        # if order == 0 sort by blocknr
        if order <= len(self.genomes) and order > -1:
            sortedLCBs = sorted(self.LCBs, key=lambda lcb: lcb.number)
                
            if order > 0:
                sortedLCBs = sorted(sortedLCBs, key=lambda lcb, order=order: 
                                    lcb.getEntry(order).start if lcb.getEntry(order) is not None else sys.maxsize
                                   )
            return(sortedLCBs)
        else:
            raise ParameterError("order", order, "between 0 and " + str(len(self.genomes)) + " (number of genomes in XMFA)")
    
    def isInvalid(self):
        for genomenr in self.genomes.keys():
            ordered = self.getSortedLCBs(genomenr)
            lastend = 0
            invalid = False
            for lcb in ordered:
                entry = lcb.getEntry(genomenr)
                if entry is None:
                    break
                    
                invalid = ( entry.start - lastend != 1 )
                if invalid: 
                    break
                    
                lastend = entry.end
                
            if invalid:
                break
        return invalid
                
                
    
        
class LCB:
    
    _iupac_dict = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "B": "CGT",
            "D": "AGT",
            "H": "ACT",
            "K": "GT",
            "M": "AC",
            "R": "AG",
            "S": "CG",
            "V": "ACG",
            "W": "AT",
            "Y": "CT",
            "AC": "M",
            "AG": "R",
            "AT": "W",
            "CG": "S",
            "CT": "Y",
            "GT": "K",
            "ACG": "V",
            "ACT": "H",
            "AGT": "D",
            "CGT": "B",
            "ACGT": "N",
            "-": "-",
            "N": "ACGT"
        }
    
    
    def __init__(self, number=0):
        self.number = number
        self.entries = []
        self.length = 0
        random.seed(a=RANDOM_SEED)
        
        
    def addEntries(self, entries):
        if type(entries) is not list:
            entries = [entries]
        
        length = len(entries[0].sequence)
        
        # check if new entries have same length as already added ones 
        # and whether all new entries have the same length
        if (self.length > 0 and length != self.length) or (not all(len(e.sequence) == length for e in entries)):
            raise LcbInputError(self.number)
       
        self.length = length
        self.entries.extend(entries)
        
        
    def getEntry(self, genomeNr):
        entry = None
        for e in self.entries:
            if e.genomeNr == genomeNr:
                entry = e
                break
        return entry

    
    def reverseComplementEntries(self):
        for e in self.entries:
            e.reverseComplement()
        
        
    def consensusSequence(self, unambiguous):
        if len(self.entries) > 0:
            if unambiguous:
                repl_method = self._random
            else:
                repl_method = self._IUPAC
                                
            consseq = ''.join([ repl_method((lambda i=i: [e.sequence[i] for e in self.entries])()) 
                                for i in range(len(self.entries[0].sequence))
                            ])  # get i'th letter in sequence of each entry of LCB   
            return consseq
        else:
            return ""
    
        
    def _IUPAC(self, bases):
        # steps:
        # convert characters into unambiguous code (A,C,G,T)
        # remove duplicates from list
        # sort list
        # remove gap character
        # convert list of bases to one ambiguous character (A,C,G,T,M,R,W,S,Y,K,V,H,D,B,N)
        # return character 
        try:
            c = ''.join(sorted(list(set(''.join([ self._iupac_dict[x.upper()] for x in bases]))))).replace("-", "")
        except KeyError as e:
            raise ConsensusFastaInputError(e.args[0])
        else:
            try:
                return self._iupac_dict[c]
            except KeyError as e:
                raise ConsensusFastaInputError(e.args[0])

    
    def _random(self, bases):
        bases = ''.join(bases).replace("-", "").upper()
        bases = bases.replace('N', "")
        if len(bases) == 0:
            return "N"
        else:
            random_choice = bases[random.randint(0,len(bases)-1)].upper()
            if not (random_choice in ('A', 'C', 'G', 'T')):
                try:
                    unambiguous_random = self._iupac_dict[random_choice]
                except KeyError as e:
                    raise ConsensusFastaInputError(e.args[0])
                else:
                    return unambiguous_random[random.randint(0, len(unambiguous_random)-1)]
            else:
                return random_choice
    
class SequenceEntry:
    
    def __init__(self, genomeNr, start, end, strand, sequence):
        self.genomeNr = int(genomeNr)
        self.sequence = sequence
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        
    
    def __setattr__(self, name, value):
        self.__dict__[name] = value
        if name == "sequence":
            self._gapDict()
    
    
    def _gapDict(self):
        self.gaps = {i.start() : i.end() for i in re.finditer('-+', self.sequence)}

        
    def getSubGapList(self, start=0, end=None):
            
        if end is None:
            end = len(self.sequence)
            
        subgaps = {gstart: gend for gstart, gend in self.gaps.items() if gend > start and gstart <= end }
        if len(subgaps) > 0:
            sortedKeys = sorted(subgaps.keys())
            smallestStart = sortedKeys[0]
            highestStart = sortedKeys[-1]
            
            # if a gap streches over borders of region, make region borders new start or end of gap
            if subgaps[highestStart] > end:
                subgaps[highestStart] = end + 1
            if smallestStart < start:
                subgaps[start] = subgaps[smallestStart]
                del subgaps[smallestStart]
        
        return subgaps

        
    def reverseComplement(self):
        self.strand = ("+" if self.strand == "-" else "-")
        seq = Seq(self.sequence)
        self.sequence = str(seq.reverse_complement())

        
    def getPositionWithinEntryWithGaps(self, posWithinBlockWithoutGaps):
        if posWithinBlockWithoutGaps == 1:
            return posWithinBlockWithoutGaps 
            
        curNrOfNonGaps = 0
        posWithinBlock = 0
        lastGapEnd = 0
    
        # go through gaps and add up non-gap characters 
        # stop if wanted position is located before next gap
        # calculate final position
        
        for start, end in sorted(self.gaps.items()):
            newNrOfNonGaps = curNrOfNonGaps + (start - posWithinBlock)
            if newNrOfNonGaps > posWithinBlockWithoutGaps:
                break
            else:
                posWithinBlock = end
                curNrOfNonGaps = newNrOfNonGaps
                
        posWithinBlock += (posWithinBlockWithoutGaps - curNrOfNonGaps)
        
        return posWithinBlock
        
    
class Consensus:
    
    def __init__(self, sequence="", order=0, xmfaFile="", fastaFile=""):
        self.sequence = sequence
        self.order = int(order)
        self.xmfaFile = os.path.abspath(xmfaFile)
        self.fastaFile = os.path.abspath(fastaFile)
     
    
    def fromAlignment(self, alignment, order, fastaFile, unambiguous):
        self.order = int(order)
        self.xmfaFile = alignment.xmfaFile
        self.fastaFile = os.path.abspath(fastaFile)
        self.sequence, self.blockStartIndices = alignment.getConsensus(unambiguous, order)
        
        
    def getFastaHeader(self, name):
        header = name + ";" + str(self.order) + "|" + self.xmfaFile
        return(header)

        
    def getUndelimitedSequence(self):
        seq = self.sequence
        parts = [seq[i:j] for i,j in zip(self.blockStartIndices, self.blockStartIndices[1:]+[None])]
        delLen = len(BLOCK_DELIMITER)
        
        seq = ''.join([part[:-delLen] for part in parts[:-1]] + [parts[-1]])
        
        return(seq)