import argparse
import os 
import sys
import traceback
import re
import itertools
import pdb
import bisect
import collections
import copy

from Bio import pairwise2
from Bio.Seq import Seq

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
        self.message = "ERROR: XMFA with more than 2 genomes provided for splitting or merging LCBs. Please align genomes to consensus sequence one by one, creating a new consensus sequence for every added genome."
        
        
class ConsensusFastaIdxFormatError(InputError):
    
    def __init__(self, message):
        self.message = "ERROR: Format of consensus fasta index file not correct:" + message

        
class CoordinateOutOfBoundsError(InputError):

    def __init__(coord, source):
        source = "consensus" if source == "c" else str(source)
        self.message = "ERROR: Position (" + str(coord) + ") not part of sequence (" + source + "."
        

class CoordinatesInputError(InputError):
    
    def __init__():
        self.message = ("ERROR: Please provide indices for mapping in correct format: "
                       "First line: source_seq\tdest_seq[,dest_seq2,...] using \"c\" or sequence number. Then one coordinate per line.")
    
        
        
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
        
    def getConsensus(self, order=0):
        sortedLCBs = self.getSortedLCBs(order)    
        delim = Parser().blockDelimiter
        
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
        
        c = ''.join(sorted(list(set(''.join([ self._iupac_dict.get(x.upper(), "|") for x in bases]))))).replace("-", "")
        return self._iupac_dict.get(c, "|")

        
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
        
    
    def fromAlignment(self, alignment, order, fastaFile):
        self.order = int(order)
        self.xmfaFile = alignment.xmfaFile
        self.fastaFile = os.path.abspath(fastaFile)
        self.sequence = alignment.getConsensus(order)
        
        
    def getFasta(self, name, delimiter=True):
        header = name + ";" + str(self.order) + "|" + self.xmfaFile
        
        seq = self.sequence
        if not delimiter:
            seq = seq.replace(Parser().blockDelimiter, '')
        
        return (">"+header+"\n"+"\n".join(re.findall(".{1,80}", seq))+"\n")   

        
    def _getDelimiterPositions(self):
        delimiter = Parser().blockDelimiter
        self.blockStartIndices = [ (i.end()) for i in re.finditer(delimiter, self.sequence)] # store end coordinate of delimiter
        self.blockStartIndices[:0] = [0]

        
class Realigner:
    
    def realign(self, alignment):
        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()
        
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
        
    
class Mapper:
    
    def mapCoordinates(self, alignment, consensus, source, dests, coordinates):
        
        if type(dests) is not list:
            dests = [dests]

        dests = sorted(dests)
        coordinates = sorted(coordinates)
        coord_dict = collections.defaultdict(dict)
        
        addC = ("c" in dests)
        
        # remove consensus-dest from list (calculation is different)
        if addC:
            dests = sorted(dests)[:-1]
        
        if source == "c":
        
            for coord in coordinates:
                
                idx = bisect.bisect_left(consensus.blockStartIndices, coord) 
                idx -= 1
                lcb = alignment.LCBs[idx]
                
                # calculate position within current consensus block - there are no gaps in consensus sequence!
                posWithinBlock = coord - consensus.blockStartIndices[idx] - 1
                
                if addC:
                    coord_dict[coord]["c"] = coord
                
                coord_dict[coord].update(self._getCoordsForEntries(lcb.entries, dests, posWithinBlock))
                            
        else:
            
            sourceBlocks = {}
            # store ends of blocks for lcb finding lcb for coords
            for i in range(len(alignment.LCBs)):
                e = alignment.LCBs[i].getEntry(int(source))
                if e is not None:
                    sourceBlocks[e.end] = {"lcb":i, "entry":e }
            
            for coord in coordinates:
                sourceEnds = sorted(sourceBlocks.keys())
                idxInEnds = bisect.bisect_left(sourceEnds, coord)
                lcbIdx = sourceBlocks[sourceEnds[idxInEnds]]["lcb"]
                
                if lcbIdx > len(alignment.LCBs):
                    raise CoordinateOutOfBoundsError(coord, source)
    
                sourceEntry = sourceBlocks[sourceEnds[idxInEnds]]["entry"]
                lcb = alignment.LCBs[lcbIdx]
                
                if sourceEntry.strand == "+":
                    posWithinBlockWithoutGaps = coord - sourceEntry.start
                else:
                    posWithinBlockWithoutGaps = sourceEntry.end - coord
            
                # add consensus coordinates to dict
                if addC:
                    consLength = sum([lcb.length for lcb in alignment.LCBs[0:lcbIdx]])
                    coord_dict[coord]["c"] = consLength + posWithinBlockWithoutGaps + 1
                    
                # check if dests other than consensus are needed
                if len(dests) > 0:
                
                    curNrOfNonGaps = 0
                    posWithinBlock = 0
                    lastGapEnd = 0
                
                    # go through gaps and add up non-gap characters 
                    # stop if wanted position is located before next gap
                    # calculate final position
                    # coordinate is already strand specific
                    
                    for start, end in sourceEntry.gaps.items():
                        newNrOfNonGaps = curNrOfNonGaps + (start - posWithinBlock)
                        if newNrOfNonGaps >= posWithinBlockWithoutGaps:
                            break
                        else:
                            posWithinBlock = end
                            curNrOfNonGaps = newNrONonGaps
                            
                    posWithinBlock += (posWithinBlockWithoutGaps - curNrOfNonGaps)
                    coord_dict[coord].update(self._getCoordsForEntries(lcb.entries, dests, posWithinBlock))
            

        return coord_dict
        
        
    def _getCoordsForEntries(self, entries, dests, posWithinBlock):
        coord_dict = {}
        
        for e in entries: # faster than looping through entries everytime to get entry of genome x
            if str(e.genomeNr) in dests: 
            
                if e.strand == "+":
                    eSubgaps = e.getSubGapList(0,posWithinBlock)
                    eSumgaps = 0
                    if len(eSubgaps) > 0:
                        eSumgaps = sum(end-start for start, end in eSubgaps.items())
                        
                    coord_dict[str(e.genomeNr)] = e.start + (posWithinBlock - eSumgaps)
                else:
                    eSubgaps = e.getSubGapList(e.end - posWithinBlock, e.end)
                    eSumgaps = 0
                    if len(eSubgaps) > 0:
                        eSumgaps = sum(end-start for start, end in eSubgaps.items())
                        
                    coord_dict[str(e.genomeNr)] = (e.end - (posWithinBlock - eSumgaps)) * -1
        
        return coord_dict

        
class Merger:

    def mergeLCBs(self, alignment, consensusGenomeNr, newGenomeNr):
        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()
            
        lcbs = alignment.getSortedLCBs(newGenomeNr)
       
        # do not create small (less bp than 10) LCBs by splitting, but append/prepend sequence         
        mergedSplitLCBs = []
        
        for i in range(len(lcbs)):
            
            lcb = lcbs[i]
            newEntry = lcb.getEntry(newGenomeNr)
            consensusEntry = lcb.getEntry(consensusGenomeNr)
            tryNextEntry = True
            
            # check if new entry is small and only created by splitting (consensus is None)
            if consensusEntry is None and newEntry is not None and len(newEntry.sequence) <= 10:
                
                nrGaps = len(newEntry.sequence)
                #sequence = newEntry.sequence
                
                # try to append to previous entry
                if i > 0:
                    lastLCB = mergedSplitLCBs[-1]
                    lastNewEntry = lastLCB.getEntry(newGenomeNr)
                    if lastNewEntry is not None:
                        tryNextEntry = False
                        
                        # check if there is a gap at the end of former entry
                        sequence = lastNewEntry.sequence
                        m = re.search("-+$", sequence)
                        if m is not None:
                            pos = m.start(0)
                            nrGaps -= (len(sequence) - pos)
                            sequence = sequence[:pos]
                            
                        
                        lastNewEntry.sequence = sequence + newEntry.sequence
                        lastNewEntry.end = newEntry.end
                    
                        lastConsensusEntry = lastLCB.getEntry(consensusGenomeNr)
                        if lastConsensusEntry is not None:
                            lastConsensusEntry.sequence += ("-"*nrGaps)
                    
                    lastLCB.length += nrGaps
                    
                # previous entry did not work, try to prepend to next entry
                if tryNextEntry and len(lcbs) > 1:
                    nextLCB = splitLCBs[i+1]
                    nextNewEntry = lastLCB.getEntry(newGenomeNr)
                    if nextNewEntry is not None:
                        tryNextEntry = False
                        
                        # check if there is a gap at the start of former entry
                        sequence = nextNewEntry.sequence
                        m = re.search("^-+", sequence)
                        if m is not None:
                            pos = m.end(0)
                            sequence = sequence[pos:]
                            nrGaps -= pos
                        
                        nextNewEntry.sequence = newEntry.sequence + sequence
                        nextNewEntry.start = newEntry.start
                    
                        nextConsensusEntry = lastLCB.getEntry(consensusGenomeNr)
                        if nextConsensusEntry is not None:
                            nextConsensusEntry.sequence = ("-"*nrGaps) + nextConsensusEntry.sequence
                        
                        nextLCB.length += nrGaps
                        
            # entry should not be merged or could neither be appended nor prepended 
            # add LCBs to alignment as it is
            if tryNextEntry:
                mergedSplitLCBs.append(lcb)

        merged = Alignment(alignment.xmfaFile)
        for nr, genome in alignment.genomes.items():
            merged.addGenome(genome, nr)
        
        for lcb in mergedSplitLCBs:
            merged.addLCB(lcb)
                
        return merged
        
        
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
                # if not: reverse complement all entries in LCB
                if consensusEntry.strand == "-":
                    lcb.reverseComplementEntries()
                    consensusEntry = lcb.getEntry(consensusGenomeNr)
            
                # check if there are delimiters in the sequence 
                trim =  (len(consensusEntry.sequence) > len(consensusEntry.sequence.strip('N')))
                split =  re.search("N[N-]{"+str(len(Parser().blockDelimiter)-2)+",}N", consensusEntry.sequence.strip('N')) is not None 
                if split or trim:
                    splitLcbs = self._splitLCB(lcb, consensusGenomeNr, newGenomeNr) # split LCB                    
                    for slcb in splitLcbs:
                        resolved.addLCB(slcb)
                else:
                    resolved.addLCB(lcb)
            else:
                resolved.addLCB(lcb)
        
        
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
                    newEntry.genomeNr = nrGenomes
                    recalculatedLCB.addEntries(newEntry)
                consensusEntry = lcb.getEntry(consensusGenomeNr)
                if consensusEntry is not None:
                    orgEntries = self._calculateCoordinates(consensusEntry, consensus, sortedOrgLCBs)
                    try:
                        recalculatedLCB.addEntries(orgEntries)
                    except LcbInputError as e:
                        e.message = e.message + " (Error occured in recalculating coordinates step.)"
                        raise e
                
                recalculated.addLCB(recalculatedLCB)
            
            return recalculated
        
        
    def _splitLCB(self, lcb, consensusGenomeNr, newGenomeNr):
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
        delimiterPositions[:] = (lambda middleStart=middleStart: [x + middleStart for x in delimiterPositions])()
        
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
        
        newEntry = SequenceEntry(entry.genomeNr, start, end, entry.strand, seq)
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
            
            newEntry = SequenceEntry(e.genomeNr, start, end, e.strand, sequence)
            
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
                    m = re.match("#Sequence(\d+)File\s+(.+)", line) # each sequence has an associated file (can be the same for more sequences -> multifasta)
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
    
    
    def parseConsensusIndex(self, filename):
        with open(filename, "r") as input:
            line = input.readline()
            m = re.match("#Fasta\t(.+)", line)
            if m is not None:
                fastaFile = m.group(1)
            else:
                raise ConsensusFastaIdxFormatError("Wrong format of Fasta header line.")
            line = input.readline()
            m = re.match("#XMFA\t(.+)", line)
            if m is not None:
                xmfaFile = m.group(1)
            else:
                raise ConsensusFastaIdxFormatError("Wrong format of XMFA header line.")
            
            alignment = Alignment(xmfaFile)
            
            lcb = LCB()
            lcbLength = 0
            lcbEndsList = [0]
            
            line = input.readline()
            while line:
                line = line.strip()
                
                if not (line == "" or line.startswith("#")): # skip rest of header - do I need genome info?
                    
                    fields = line.split("\t")
                    id = fields[0]
                    nr = int(fields[1])
                    start = int(fields[2])
                    end = int(fields[3])
                    strand = fields[4]
                    if id == "b":
                        if lcbLength > 0:
                            lcb.length = lcbLength
                            alignment.addLCB(lcb)
                            lcb = LCB()
                        lcbLength = (end - start) + 1
                        lcbEndsList.append(end)
                    elif id == "s":
                        e = SequenceEntry(nr, start, end, strand, '')
                        
                        if len(fields) == 6:
                            gaps = fields[5]
                            gapDict = {}
                            for interval in gaps.split(";"):
                                start, end = interval.split("-")
                                gapDict[int(start)] =  int(end)
                            e.gaps = gapDict
                        lcb.entries.append(e)
                    else:
                        raise ConsensusFastaIdxFormatError("Lines can only start with 'b' or 's'.")
                    
                line = input.readline()
                
            consensus = Consensus( sequence="", order=0, xmfaFile=xmfaFile, fastaFile=fastaFile)
            consensus.blockStartIndices = lcbEndsList
                
        return (alignment, consensus)
        
        
    def parseMappingCoordinates(self, coord_f):
        with open(coord_f) as input:
            
            header = input.readline().strip()
            source, dest = header.split("\t")
            dests = dest.split(",")
            coords = [int(line.strip()) for line in input]
            
            if source == "" or dests == "" or len(dests) == 0 or len(coords) == 0:
                raise CoordinatesInputError()
        
        return source, dests, coords
        
        
class Writer:
        _mauveFormatString = "#FormatVersion Mauve1\n"
        _mauveGenomeFile = '#Sequence{0}File\t{1}\n'
        _mauveGenomeEntry = '#Sequence{0}Entry\t{1}\n'
        _mauveGenomeFormat = '#Sequence{0}Format\t{1}\n'
        _mauveBlockHeader = '> {0}:{1}-{2} {3} {4}\n'
        
        _consensusIndexFasta = '#Fasta\t{0}\n'
        _consensusIndexXmfa = '#XMFA\t{0}\n'
        _consensusIndexBlockLine = '\nb\t{0}\t{1}\t{2}\t+\n'
        _consensusIndexSequenceLine = 's\t{0}\t{1}\t{2}\t{3}\t{4}\n'
        
        _mafFormatString = "##maf version=1\n"
        _mafSequenceHeader = "\na label={0}\n"
        _mafEntryHeader = "s {0}\t{1} {2} {3} {4} {5}\n"
        
        
        def writeXMFA(self, alignment, path, name, order=0):
            with open(path+"/"+name+".xmfa", "w") as output:
                output.write(self._mauveFormatString)
                               
                for nr, genome in sorted(alignment.genomes.items()):
                    output.write(self._mauveGenomeFile.format(nr, genome.filepath))
                    if genome.entry > 0 :
                        output.write(self._mauveGenomeEntry.format(nr, genome.entry))
                    output.write(self._mauveGenomeFormat.format(nr, genome.format))
                
                sortedLCBs = alignment.getSortedLCBs(order)
                count = 0
                for lcb in sortedLCBs:
                    count += 1
                    for entry in sorted(lcb.entries, key=lambda e: e.genomeNr):
                        output.write(self._mauveBlockHeader.format( entry.genomeNr, 
                                                                   entry.start, 
                                                                   entry.end, 
                                                                   entry.strand, 
                                                                   alignment.genomes[entry.genomeNr].filepath
                                                                  )
                                    )
                        output.write("\n".join(re.findall(".{1,80}", entry.sequence))+"\n")
                    output.write("=\n")
        
        
        def writeMAF(self, alignment, path, name, order=0):
            
            with open(path+"/"+name+".maf", "w") as output:
                output.write(self._mafFormatString)
                               
                for nr, genome in sorted(alignment.genomes.items()):
                    output.write(self._mauveGenomeFile.format(nr, genome.filepath))
                    if genome.entry > 0 :
                        output.write(self._mauveGenomeEntry.format(nr, genome.entry))
                    output.write(self._mauveGenomeFormat.format(nr, genome.format))
                
                genome_lengths = {}
                sortedLCBs = alignment.getSortedLCBs(order)
                
                for lcb in sortedLCBs:
                    for entry in lcb.entries:
                        genome_lengths[entry.genomeNr] = genome_lengths.get(entry.genomeNr, 0) + (len(entry.sequence) - entry.sequence.count('-'))
                
                count = 0
                for lcb in sortedLCBs:
                    count += 1
                    output.write(self._mafSequenceHeader.format(count))
                    for entry in sorted(lcb.entries, key=lambda e: e.genomeNr):
                        output.write(self._mafEntryHeader.format(entry.genomeNr, entry.start, (entry.end - entry.start), genome_lengths[entry.genomeNr], entry.strand, entry.sequence))
                    
                
        
        def writeMappingCoordinates(self, source, dests, coords_dict, path, name):
            with open(os.path.abspath(path + "/"+ name + ".txt"), "w") as output:
                output.write(''.join([str(source)," (source)", "\t"]))
                output.write('\t'.join(dests))
                output.write("\n")
                
                for coord, cur_dict in coords_dict.items():
                    output.write(str(coord) + "\t")
                    new_coords = [str(cur_dict[dest]) for dest in dests]
                    output.write('\t'.join(new_coords))
                    output.write("\n")

        
        def writeConsensus(self, alignment, path, name, order=0):
            filename = os.path.abspath(path+"/"+name+"_consensus.fasta")
            
            consensus = Consensus()
            self._writeConsensusIndex(alignment, filename, order)
            consensus.fromAlignment(alignment, order, filename)
            
            with open(filename, "w") as output:
                output.write(consensus.getFasta(name, False))
            
            with open(filename + ".blockseparated.fasta", "w") as output:
                output.write(consensus.getFasta(name, True))

                
        def _writeConsensusIndex(self, alignment, fastafile, order=0):
            with open(fastafile+".idx", "w") as output:
                output.write(self._consensusIndexFasta.format(fastafile))
                output.write(self._consensusIndexXmfa.format(alignment.xmfaFile))
                
                for nr, genome in sorted(alignment.genomes.items()):
                    output.write(self._mauveGenomeFile.format(nr, genome.filepath))
                    if genome.entry > 0 :
                        output.write(self._mauveGenomeEntry.format(nr, genome.entry))
                    output.write(self._mauveGenomeFormat.format(nr, genome.format))
                
                sortedLCBs = alignment.getSortedLCBs(order)
            
                consensusEnd = 0
                counter = 1
            
                for lcb in sortedLCBs:
                    output.write(self._consensusIndexBlockLine.format(counter, consensusEnd+1, consensusEnd+lcb.length))
                    consensusEnd += (lcb.length)
                    
                    counter += 1
                    for entry in lcb.entries:
                        output.write(self._consensusIndexSequenceLine.format
                                            ( entry.genomeNr, entry.start, entry.end, entry.strand,
                                              ';'.join(['-'.join( [str(start+1), str(end)]) for start, end in sorted(entry.gaps.items()) ])
                                            )
                                    )
                                            
def main():
    global args
    
    parser = Parser()
    writer = Writer()
    resolver = Resolver()
    realigner = Realigner()
    mapper = Mapper()
    merger = Merger()
    
    if args.task == "map":
        try:
            sparse_align, sparse_consensus = parser.parseConsensusIndex(args.consensus_f+".idx")
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
                    realign = realigner.realign(align)
                    writer.writeXMFA(realign, args.output_p, args.output_name + "_realign", args.order)
                
                if args.task == "merge":
                    merged = merger.mergeLCBs(align, 1, 2)
                    writer.writeXMFA(merged, args.output_p, args.output_name + "_merge", args.order)
                    
                if args.task == "consensus":
                    writer.writeConsensus(align, args.output_p, args.output_name, args.order)
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
                
                elif args.task == "maf":
                    writer.writeMAF(align, args.output_p, args.output_name, args.order)
                    
            except ParameterError as e:
                print('ERROR: Problem with parameter "{0}": Value should be {1}, but was "{2}".'.format(e.parameter, e.rangetext, e.value))
    return(0)
        
if __name__ == '__main__':
    #try:
                
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xmfa", dest="xmfa_f", help="XMFA input file")
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name", help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    parser.add_argument("-c", "--consensus", dest="consensus_f", help="consensus FASTA file used in XMFA", required=False)
    parser.add_argument("-o", "--order", dest="order", type=int, default=0, help="ordering of output (0,1,2,...) [default: %(default)s]", required=False)
    parser.add_argument("-t", "--task", dest="task", default="consensus", help="what to do (consensus|split|realign|xmfa|map|merge|maf) [default: %(default)s]", choices=["consensus", "split", "realign", "xmfa", "maf", "map", "merge"], required=False)
    parser.add_argument("-i", "--index", dest="coord_f", help="file with indices to map. First line: source_seq\tdest_seq[,dest_seq2,...] using \"c\" or sequence number. Then one coordinate per line.")
    
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
            parser.error("the following arguments are required: -x/--xmfa")
    
    main()
    #    sys.exit(0)
        
  #  except KeyboardInterrupt as e: # Ctrl-C
  #      raise e
  #  except SystemExit as e: # sys.exit()
  #      raise e
  #  except Exception as e:
  #      print('ERROR, UNEXPECTED EXCEPTION')
  #      print(str(e))
  #      traceback.print_exc()
  #      os._exit(1)