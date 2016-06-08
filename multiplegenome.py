import argparse
import os 
import sys
import traceback
import re
import itertools
import pdb
import bisect
import collections

from Bio import pairwise2
from Bio.Seq import Seq
#sys.path.append('/home/jandrasitsc/software/bin/python/lib/python3.4/')

#import swalign

class FormatError(Exception):
    pass

class InputError(Exception):
    pass

class ParameterError(Exception):
    
    def __init__(self, parameter, value, rangetext):
        self.parameter = parameter
        self.value = value
        self.rangetext = rangetext
    
class ConsensusFastaFormatError(FormatError):
    def __init__(self):
        self.message = "ERROR: Wrong format of consensus fasta header. Please rebuild consensus fasta with current script version."
       

class XMFAHeaderFormatError(FormatError):
    def __init__(self, header):
        self.message = "ERROR: XMFA sequence header \"" + header + "\" does not apply to XMFA header format rule: \"> seq:start-end strand\"."
       
       
class LcbInputError(InputError):
    
    def __init__(self, lcbnr):
        self.message = 'ERROR: Problem with LCB Nr. {0}: Entries must not be of different length.'.format(lcbnr)
        
        
class ConsensusXMFAInputError(InputError):
    
    def __init__(self):
        self.message = "ERROR: XMFA with more than 2 genomes provided for splitting LCBs. Please align genomes to consensus sequence one by one, creating a new consensus sequence for every added genome."
        
        
class Genome:

    def __init__(self, filepath, format, entry=-1):
        self.filepath = os.path.abspath(filepath)
        self.format = format
        self.entry = int(entry)
                
class Alignment:

    def __init__(self, xmfaFile):
        self.xmfaFile = os.path.abspath(xmfaFile)
        self.LCBs = []
        self.genomes = {}
        
    def addGenome(self, genome, number):
        self.genomes[int(number)] = genome
        
    def addLCBentries(self, entries):
        number = len(self.LCBs) + 1
        lcb = LCB(number)
        lcb.addEntries(entries)
        self.LCBs.append(lcb)
        
    def addLCB(self, lcb):
        number = len(self.LCBs) + 1
        lcb.number = number
        self.LCBs.append(lcb)
        
    def getConsensus(self, order=0, noDelimiter=False):
        
        sortedLCBs = self.getSortedLCBs(order)    
        delim = Parser().blockDelimiter
        
        if noDelimiter:
            delim = ""
        
        consseq = delim.join( lcb.consensusSequence() for lcb in sortedLCBs )
        
        return consseq
            

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
        
    def addEntries(self, entries):
        if type(entries) is not list:
            entries = [entries]
        
        length = len(entries[0].sequence)
        
        # check if new entries have same length as already added 
        # and whether all new entries have the same length
        if (self.length > 0 and length != self.length) or (not all(len(e.sequence) == length for e in entries)):
            raise LcbInputError(self.number)
       
        self.length = length
        self.entries.extend(entries)
        
    def getEntry(self, sequenceNr):
        entry = None
        for e in self.entries:
            if e.sequenceNr == sequenceNr:
                entry = e
                break
        return entry

    
    def reverseComplementEntries(self):
        for e in self.entries:
            e.reverseComplement()
        
        
    def consensusSequence(self):
        if len(self.entries) > 0:
            consseq = ''.join([ self._IUPAC((lambda i=i: [e.sequence[i] for e in self.entries])()) 
                                for i in range(len(self.entries[0].sequence))
                              ])   # get i'th letter in sequence of each entry of LCB
            return consseq
        else:
            return ""
            
    def _IUPAC(self, bases):
        # steps:
        # convert character into unambiguous code (A,C,G,T)
        # remove duplicates from list
        # sort list
        # remove gap character
        # convert list of bases to one ambiguous character (A,C,G,T,M,R,W,S,Y,K,V,H,D,B,N)
        # return character (or "|" if something goes wrong
        c = ''.join(sorted(list(set(''.join([self._iupac_dict.get(x, "|") for x in bases]))))).replace("-", "")
        return self._iupac_dict.get(c, "|")

        
class SequenceEntry:
    
    def __init__(self, sequenceNr, start, end, strand, sequence):
        self.sequenceNr = int(sequenceNr)
        self.sequence = sequence
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        
    
    def __setattr__(self, name, value):
        self.__dict__[name] = value
        if name == "sequence":
            self._gapList()
    
    def _gapList(self):
        self.gaps = {i.start() : i.end() for i in re.finditer('-+', self.sequence)}

    def getSubGapList(self, start=0, end=None):
            
        if end is None:
            end = len(self.sequence)
            
        subgaps = {gstart: gend for gstart, gend in self.gaps.items() if gend > start and gstart < end }
        if len(subgaps) > 0:
            sortedKeys = sorted(subgaps.keys())
            smallestStart = sortedKeys[0]
            highestStart = sortedKeys[-1]
            
            # if a gap streches over borders of region, make region borders new start or end of gap
            if subgaps[highestStart] > end:
                subgaps[highestStart] = end
            if smallestStart < start:
                subgaps[start] = subgaps[smallestStart]
                del subgaps[smallestStart]
        
        return subgaps
        
    def reverseComplement(self):
        self.strand = ("+" if self.strand == "-" else "-")
        seq = Seq(self.sequence)
        self.sequence = str(seq.reverse_complement())

        
class Consensus:
    
    def __init__(self, sequence="", order=0, xmfaFile="", fastaFile=""):
        self.sequence = sequence
        self.order = int(order)
        self.xmfaFile = os.path.abspath(xmfaFile)
        self.fastaFile = os.path.abspath(fastaFile)
        
    
    def __setattr__(self, name, value):
        self.__dict__[name] = value
        if name == "sequence":
            self._getDelimiterPositions()
        
    
    def fromAlignment(self, alignment, order, fastaFile, nodelimiter=False):
        self.order = int(order)
        self.xmfaFile = alignment.xmfaFile
        self.fastaFile = os.path.abspath(fastaFile)
        self.sequence = alignment.getConsensus(order, nodelimiter)
        
    def getFasta(self, name):
        header = name + ";" + str(self.order) + "|" + self.xmfaFile
        return (">"+header+"\n"+"\n".join(re.findall(".{1,80}", self.sequence))+"\n")   

    def _getDelimiterPositions(self):
        delimiter = Parser().blockDelimiter
        self.blockStartIndices = [ (i.end()) for i in re.finditer(delimiter, self.sequence)] # store end coordinate of delimiter
        self.blockStartIndices[:0] = [0]

        
class Realigner:
    
    #def __init__(self):
    
    def realign(self, alignment):
        # start with 2 sequences -> DO something for more or add exception (as with resolve)
        
        realigned = Alignment(alignment.xmfaFile)
        for nr, genome in alignment.genomes.items():
            realigned.addGenome(genome, nr)
        
        
        # go through lcbs, skip one-entry ones
        for lcb in alignment.getSortedLCBs(0):
            if len(lcb.entries) == 1:
                realigned.addLCB(lcb)
            else:
                entryOne = lcb.entries[0]
                entryTwo = lcb.entries[1]
                
                # get regions to realign
                oneFirstTwoSecond = self._getRealignRegions(entryOne.getSubGapList(), entryTwo.getSubGapList()  )
                
                if len(oneFirstTwoSecond) > 0:
                    seqOne, seqTwo = self._realign(entryOne.sequence, entryTwo.sequence, oneFirstTwoSecond)
                    entryOne.sequence = seqOne
                    entryTwo.sequence = seqTwo
                
                # get regions to realign for updated entries with second entry first
                twoFirstOneSecond = self._getRealignRegions(entryTwo.getSubGapList(), entryOne.getSubGapList() )
                
                if len(twoFirstOneSecond) > 0:
                    seqTwo, seqOne = self._realign(entryTwo.sequence, entryOne.sequence, twoFirstOneSecond)
                    entryOne.sequence = seqOne
                    entryTwo.sequence = seqTwo
                
                newlcb = LCB()
                newlcb.addEntries([entryOne, entryTwo])
                
                realigned.addLCB(newlcb)
                
        return realigned
    
    
    def _getRealignRegions(self, gapsForStart, gapsForEnd):
        endsDict = { end : start for start, end in gapsForEnd.items() }
        locList = [start for start, end in gapsForStart.items()] + list(endsDict.keys())
        
        regionStarts = [item for item, count in collections.Counter(locList).items() if count > 1]
        
        regions = []
        for start in regionStarts:
            regions.append( [(start, gapsForStart[start]), (endsDict[start], start)])
        
        return regions

    
    def _realign(self, seqOne, seqTwo, realignRegions):
    
        realignRegions = sorted(realignRegions)
    
        index_offset = 0
    
        for interval in realignRegions:
            minIndex, minSeqLength = min(enumerate( [interval[0][1] - interval[0][0], interval[1][1] - interval[1][0] ] ), key=lambda p: p[1])
            
            if minSeqLength < 10:
            
                seqStart = interval[minIndex][0] - index_offset - minSeqLength
                seqEnd = interval[minIndex][1] -index_offset + minSeqLength
                
                alignments = pairwise2.align.globalxx(seqOne[seqStart:seqEnd].replace("-", "") , seqTwo[seqStart:seqEnd].replace("-", ""))
                
                
                maxscore = max( [x[2] for x in alignments] )
                alignments = (lambda maxscore=maxscore: [item for item in alignments if item[2] == maxscore])()
                
                minlength = min( [x[4] for x in alignments] )
                alignments = (lambda minlength=minlength: [item for item in alignments if item[4] == minlength])()
                
                seqOne = alignments[0][0].join([ seqOne[:seqStart], seqOne[seqEnd:] ])
                seqTwo = alignments[0][1].join([ seqTwo[:seqStart], seqTwo[seqEnd:] ])
                
                index_offset += ((seqEnd - seqStart) - minlength)

        
        return (seqOne, seqTwo)
        
class Resolver:        
    
    def resolveMultiAlignment(self, alignment, consensus, orgAlignment):
        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()
    
        consensusGenomeNr = ( 1 if alignment.genomes[1].filepath == consensus.fastaFile else 2)
        newGenomeNr = (1 if consensusGenomeNr == 2 else 2)
    
        resolved = Alignment(alignment.xmfaFile)
        for nr, genome in alignment.genomes.items():
            resolved.addGenome(genome, nr)
        
        for lcb in alignment.LCBs:
            consensusEntry = lcb.getEntry(consensusGenomeNr)
            if consensusEntry is not None:
                # consensus sequence entry of LCB should be on forward strand for easier calculation of coordinates
                # if not reverse complement all entries in LCB
                if consensusEntry.strand == "-":
                    lcb.reverseComplementEntries()
                    consensusEntry = lcb.getEntry(consensusGenomeNr)
            
                # check if there are delimiters in the sequence 
                trim =  (len(consensusEntry.sequence) > len(consensusEntry.sequence.strip('N')))
                split =  re.search("N[N-]{"+str(len(Parser().blockDelimiter)-2)+",}N", consensusEntry.sequence.strip('N')) is not None 
                if split or trim:
                    splitLcbs = self._splitLCB(lcb, consensusGenomeNr) # split LCB
                    for slcb in splitLcbs:
                        resolved.addLCB(slcb)
                else:
                    resolved.addLCB(lcb)
            else:
                resolved.addLCB(lcb)
        
        #Writer().writeXMFA(resolved, "/home/jandrasitsc/TB_analysis/genomes/mauve/parsertest/testsequence", "unrecalculated", 0)
        
        recalculated = Alignment(orgAlignment.xmfaFile)
        for nr, genome in orgAlignment.genomes.items():
            recalculated.addGenome(genome, nr)
        nrGenomes = len(orgAlignment.genomes)+1
        recalculated.addGenome(alignment.genomes[newGenomeNr], nrGenomes)
        try:
            sortedOrgLCBs = orgAlignment.getSortedLCBs(consensus.order)    
        except ParameterError:
            raise ConsensusFastaFormatError()
        else:
            for lcb in resolved.LCBs:
                recalculatedLCB = LCB(lcb.number)
                newEntry = lcb.getEntry(newGenomeNr)
                
                if newEntry is not None:
                    newEntry.sequenceNr = nrGenomes
                    recalculatedLCB.addEntries(newEntry)
                consensusEntry = lcb.getEntry(consensusGenomeNr)
                if consensusEntry is not None:
                    #pdb.set_trace()
                    orgEntries = self._calculateCoordinates(consensusEntry, consensus, sortedOrgLCBs)
                    try:
                        recalculatedLCB.addEntries(orgEntries)
                    except LcbInputError as e:
                        print(e.message + " (Error occured in recalculating coordinates step.)")
                
                recalculated.addLCB(recalculatedLCB)
            
            return recalculated
        
        
    def _splitLCB(self, lcb, consensusGenomeNr):
        splitLCBs = []
        delimiterPositions = []
        startPositions = [0]
        
        consensusEntry = lcb.getEntry(consensusGenomeNr)
        
        middleStart = 0
        middleEnd = lcb.length
        
        # check if there is a delimiter sequence at the start
        # and search from end of match afterwards
        m = re.search("^N[N-]*", consensusEntry.sequence)
        if m is not None:
            startPositions.append(m.end(0))
            middleStart = m.end(0)
        
        # check if there is a delimiter sequence at the end
        # and search until start of match afterwards
        m = re.search("[N-]*N$", consensusEntry.sequence[middleStart:])
        if m is not None:
            delimiterPositions.append([m.start(0)])
            middleEnd = m.start(0) + middleStart
        
        
        # N stretches in the middle must have delimiter length and be between start and end delimiter sequence
        delimiterPositions.extend([ (m.start(0), m.end(0)) 
                                    for m in re.finditer("N[N-]{"+str(len(Parser().blockDelimiter)-2)+",}N", 
                                                         consensusEntry.sequence[middleStart:middleEnd]
                                                        )
                                  ] ) 
        
        delimiterPositions = sorted(list(itertools.chain.from_iterable(delimiterPositions)))
        delimiterPositions[:] = (lambda start=middleStart: [x + start for x in delimiterPositions])()
        
        if len(delimiterPositions) == 0 and len(startPositions) < 2: # no entries with delimiter sequence in LCB
            return [lcb]
        
        indices = startPositions + delimiterPositions +[lcb.length]
        
        for i in range(len(indices)-1):
            entries = list(filter(None, (lambda i=i, lcb=lcb, self=self, indices=indices: 
                                            [ self._getNewEntry(e, indices[i], indices[i+1]) 
                                              for e in lcb.entries
                                            ]
                                        )()
                                  ) # add all entries that are not None
                          )
            if len(entries) > 0:
                if len(entries) == 1: # only one entry - sequence of second genome not present in current LCB --> remove all gaps
                    entries[0].sequence = entries[0].sequence.replace("-", "") 
                slcb = LCB()
                try:
                    slcb.addEntries(entries)
                except LcbInputError as e:
                    raise LcbInputError(e.message + " (Error occured in splitting step.)")
                else:    
                    splitLCBs.append(slcb)
        
        return splitLCBs
        

    def _getNewEntry(self, entry, splitstart, splitend):
        seq = entry.sequence[splitstart:splitend]
        if seq == "" or re.search("[^N-]", seq) is None:
            return None
        splitlen = splitend - splitstart
        
        sumgaps = 0
        
        subgaps = entry.getSubGapList(0, splitstart)
        if len(subgaps) > 0:
            sumgaps = sum(end-start for start, end in subgaps.items())
        
        start = entry.start + splitstart - sumgaps
        end = start + splitlen - 1 
        
        newEntry = SequenceEntry(entry.sequenceNr, start, end, entry.strand, seq)
        newEntry.end = newEntry.end - sum( end-start for start, end in newEntry.gaps.items() )
        
        return newEntry

    
    
    def _calculateCoordinates(self, consensusEntry, consensus, orgLCBlist):
        
        # calculate in which consensus block this entry is located
        idx = bisect.bisect_left(consensus.blockStartIndices, consensusEntry.start) 
        
        idx -= 1
        orgLCB = orgLCBlist[idx]
        
        # calculate start and end of sequence within current consensus block
        startWithinBlock = consensusEntry.start - consensus.blockStartIndices[idx] - 1
        endWithinBlock = consensusEntry.end - consensus.blockStartIndices[idx] - 1
        
        # count gaps in consensus entry
        sumgaps = 0
        if len(consensusEntry.gaps) > 0:
            sumgaps = sum(end-start for start, end in consensusEntry.gaps.items())
        
        # get sequence of original blocks the same length as the consensus entry sequence minus number of gaps
        endSequenceSub = startWithinBlock + len(consensusEntry.sequence) - sumgaps
        
        orgEntries = []
        for e in orgLCB.entries:
            start = e.start + startWithinBlock
            end = e.start + endWithinBlock
            sequence = e.sequence[startWithinBlock:endSequenceSub]
            
            # count gaps in original sequence before start of sequence to get correct start and end for entry
            eSumgaps = 0
            eSubgaps = e.getSubGapList(0,startWithinBlock)
            if len(eSubgaps) > 0:
                eSumgaps = sum(end-start for start, end in eSubgaps.items())
            
            start = start - eSumgaps
            end = end - eSumgaps
            
            newEntry = SequenceEntry(e.sequenceNr, start, end, e.strand, sequence)
            
            newEntry.end = newEntry.end - sum( end-start for start, end in newEntry.gaps.items() )
            
            # include all gaps in consensus sequence in original sequence for correct alignment
            for gstart, gend in sorted(consensusEntry.gaps.items()):
                sequence = self._insertGap(sequence, gstart, gend-gstart)
            
            if sequence != "" and re.search("[^N-]", sequence) is not None:
                newEntry.sequence = sequence
                orgEntries.append(newEntry)
        
        return orgEntries
    
    
    def _insertGap(self, sequence, position, length):
        seq = ("-"*length).join([sequence[:position], sequence[position:]])
        return seq

        
class Parser:
    blockDelimiter = 'N' * 1000

    def parseXMFA(self, filename):
        alignment = Alignment(filename)
        
        with open(filename, "r") as xmfa:
            line = xmfa.readline()
            seq = ""
            start = 0
            end = 0
            seqNr = 0
            ses = []
            while line:
                line = line.rstrip()
                if line.startswith("#"): #parse comment section
                    m = re.match("#Sequence(\d)File\s+(.+)", line) # each sequence has an associated file (can be the same for more sequences -> multifasta)
                    if m is not None:
                        number = m.group(1)
                        fn = m.group(2)
                        entry = -1
                        line = xmfa.readline()
                        m = re.match("#Sequence"+number+"Entry\s+(\d+)", line) # with multifasta files sequence - entry numbers are reported in line after filename
                        if m is not None:
                            entry = m.group(1)
                            line = xmfa.readline()
                            
                        m = re.match("#Sequence"+number+"Format\s+(\w+)", line) 
                        if m is not None:
                            format = m.group(1)
                            line = xmfa.readline()
                        genome = Genome(fn, format, entry)
                        alignment.addGenome(genome, number)
                        continue
                
                elif line.startswith(">"): # a sequence start was encountered
                    if len(seq) > 0: # save previous sequence
                        ses.append(SequenceEntry(seqNr, start, end, strand, seq))
                        seq = ""
                    m = re.match(">\s*(\d+):(\d+)-(\d+) ([+-]) ", line)
                    if m is not None:
                        seqNr = m.group(1)
                        start = m.group(2)
                        end = m.group(3)
                        strand = m.group(4)
                    else:
                        raise XMFAHeaderFormatError(line.strip())
                elif line.startswith("="):
                    ses.append(SequenceEntry(seqNr, start, end, strand, seq))
                    
                    alignment.addLCBentries(ses)
                    
                    seq = ""
                    ses = []
                else:
                    seq += line
                    
                line = xmfa.readline()
                    
        return alignment
    
    
    def parseConsensus(self, filename):
        
        with open(filename, "r") as input:
            line = input.readline()
            m = re.match("^>[^;]+;(\d+)\|(.*)", line)
            if m is not None:
                order = m.group(1)
                xmfaFile = m.group(2)
            else:
                raise ConsensusFastaFormatError()
            line = input.readline()
            sequence = ""
            while line:
                sequence = sequence + line.strip()
                line = input.readline()
            
            try:
                cons = Consensus(sequence, order, xmfaFile, filename)
            except ParameterError:
                raise ConsensusFastaFormatError()
                
            return cons
    
    

class Writer:
        mauveFormatString = "#FormatVersion Mauve1\n"
        mauveGenomeFile = '#Sequence{0}File\t{1}\n'
        mauveGenomeEntry = '#Sequence{0}Entry\t{1}\n'
        mauveGenomeFormat = '#Sequence{0}Format\t{1}\n'
        mauveBlockHeader = '> {0}:{1}-{2} {3} {4}\n'
        
        def writeXMFA(self, alignment, path, name, order=0):
            with open(path+"/"+name+".xmfa", "w") as output:
                output.write(self.mauveFormatString)
                               
                for nr, genome in sorted(alignment.genomes.items()):
                    output.write(self.mauveGenomeFile.format(nr, genome.filepath))
                    if genome.entry > 0 :
                        output.write(self.mauveGenomeEntry.format(nr, genome.entry))
                    output.write(self.mauveGenomeFormat.format(nr, genome.format))
                
                sortedLCBs = alignment.getSortedLCBs(order)
                count = 0
                for lcb in sortedLCBs:
                    count += 1
                    for entry in sorted(lcb.entries, key=lambda e: e.sequenceNr):
                        output.write(self.mauveBlockHeader.format( entry.sequenceNr, 
                                                                   entry.start, 
                                                                   entry.end, 
                                                                   entry.strand, 
                                                                   alignment.genomes[entry.sequenceNr].filepath
                                                                  )
                                    )
                        output.write("\n".join(re.findall(".{1,80}", entry.sequence))+"\n")
                    output.write("=\n")
            
        
        def writeConsensus(self, alignment, path, name, order=0, nodelimiter=False):
            filename = path+"/"+name+"_consensus.fasta"
            consensus = Consensus()
            consensus.fromAlignment(alignment, order, filename, nodelimiter)
            with open(filename, "w") as output:
                output.write(consensus.getFasta(name))
        
def main():
    global args
    
    parser = Parser()
    writer = Writer()
    resolver = Resolver()
    realigner = Realigner()
    
    try:
        align = parser.parseXMFA(args.xmfa_f)
            
    except (XMFAHeaderFormatError, LcbInputError) as e:
        print(e.message + "(" + args.xmfa_f + ")")
    else:
        try:
            if args.task == "realign":
                realign = realigner.realign(align)
                writer.writeXMFA(realign, args.output_p, args.output_name + "_realign", args.order)
                
            if args.task == "consensus":
         #       pdb.set_trace()
                writer.writeConsensus(align, args.output_p, args.output_name, args.order, args.nodelimiter)
            elif args.task == "split":
                
                try:
                    consensus = parser.parseConsensus(args.consensus_f)
                    org_align = parser.parseXMFA(consensus.xmfaFile)
                    splitblocks_align = resolver.resolveMultiAlignment(align, consensus, org_align)
                except (XMFAHeaderFormatError, LcbInputError) as e:
                        print(e.message + "(" + consensus.xmfaFile + ")")
                except (ConsensusFastaFormatError, ConsensusXMFAInputError) as e:
                    print(e.message)
                else:
                    writer.writeXMFA(splitblocks_align, args.output_p, args.output_name+"_split", args.order)
            
            elif args.task == "xmfa":
                writer.writeXMFA(align, args.output_p, args.output_name, args.order)
                
                
        except ParameterError as e:
            print('ERROR: Problem with parameter "{0}": Value should be {1}, but was "{2}".'.format(e.parameter, e.rangetext, e.value))
        
        
if __name__ == '__main__':
    try:
                
        parser = argparse.ArgumentParser()
        parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file", required=True)
        parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
        parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
        parser.add_argument("-c", "--consensus", dest="consensus_f", help="consensus FASTA file used in XMFA", required=False)
        parser.add_argument("-o", "--order", dest="order", type=int, default=0, help="ordering of output (0,1,2,...) [default: %(default)s]", required=False)
        parser.add_argument("-t", "--task", dest="task", default="consensus", help="what to do (consensus|split|realign|xmfa) [default: %(default)s]", choices=["consensus", "split", "realign", "xmfa"], required=False)
        parser.add_argument( "--nodelimiter", dest="nodelimiter", action="store_true", help="print consensus sequence without block delimiter")
        
        args = parser.parse_args()
        
        if args.task == "split" and args.consensus_f is None:
             parser.error("Please provide a consensus-sequence file (-c/--consensus) for the \"split\"-task (-t/--task).")
            
        if args.task != "consensus" and args.nodelimiter:
            print("warnings: flag '--nodelimiter' has no effect if task is unequal to 'consensus'. " , file=sys.stderr)
        
        
        #if len(args) < 1:
        #    parser.error ('missing argument')
        #else:
        main()
        sys.exit(0)
        
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)