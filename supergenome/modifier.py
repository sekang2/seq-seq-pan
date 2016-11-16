import re
import collections

from Bio import pairwise2


from supergenome.exception import ConsensusXMFAInputError
from supergenome.base import *

class Separator:
    
    def separateLCBs(self, alignment, length):
        separated = Alignment(alignment.xmfaFile)
        for nr, genome in alignment.genomes.items():
            separated.addGenome(genome, nr)
        
        for lcb in alignment.LCBs:
            if lcb.length <= length:
                for entry in lcb.entries:
                    seq = entry.sequence.replace("-", "") 
                    newEntry = SequenceEntry(entry.genomeNr, entry.start, entry.end, entry.strand, seq)
                    separated.addLCBEntries(newEntry)
            else:
                separated.addLCB(lcb)
                
        return separated

class Merger:

    def mergeLCBs(self, alignment, consensusGenomeNr, newGenomeNr, blockLength):
        
        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()
            
        lcbs = alignment.getSortedLCBs(newGenomeNr)
       
        # do not create small (less bp than blocklength) LCBs by splitting, but append/prepend sequence         
        mergedSplitLCBs = []
        
        for i in range(len(lcbs)):
            
            lcb = lcbs[i]
            newEntry = lcb.getEntry(newGenomeNr)
            consensusEntry = lcb.getEntry(consensusGenomeNr)
            prepend = False
            append = False
            revcompl = False
            
            neighbourNewEntry = None
            neighbourConsensusEntry = None
            
            # check if new entry is small and only created by splitting or aligning of new genome (consensus is None)
            if consensusEntry is None and newEntry is not None and len(newEntry.sequence) <= blockLength:
                
                nrGaps = len(newEntry.sequence)
                
                # try to append to previous entry
                if i > 0:
                    neighbourLCB = mergedSplitLCBs[-1]
                    neighbourNewEntry = neighbourLCB.getEntry(newGenomeNr)
                    if neighbourNewEntry is not None and (newEntry.start - neighbourNewEntry.end) == 1:
                        
                        if newEntry.strand == neighbourNewEntry.strand:
                            append = True
                        else:
                            prepend = True
                            revcompl = True
                        
                # previous entry did not work, try to prepend to next entry
                if not(prepend) and not(append) and len(lcbs) > i:
                    neighbourLCB = lcbs[i+1]
                    neighbourNewEntry = neighbourLCB.getEntry(newGenomeNr)
                    if neighbourNewEntry is not None and (neighbourNewEntry.start - newEntry.end) == 1:
                        
                        if newEntry.strand == neighbourNewEntry.strand:
                            prepend = True
                        else:
                            append = True
                            revcompl = True
                        
                        
            if append or prepend:
            
                if revcompl:
                    newEntry.reverseComplement()
            
                sequence = neighbourNewEntry.sequence
                neighbourConsensusEntry = neighbourLCB.getEntry(consensusGenomeNr)
                neighbourLCB.length += nrGaps
            
                neighbourNewEntry.end = max(newEntry.end, neighbourNewEntry.end)
                neighbourNewEntry.start = min(newEntry.start, neighbourNewEntry.start)
            
                if append:
                    
                    neighbourNewEntry.sequence = sequence + newEntry.sequence
                    
                    if neighbourConsensusEntry is not None:
                        neighbourConsensusEntry.sequence += ("-"*nrGaps)
                
                elif prepend:
                    
                    neighbourNewEntry.sequence = newEntry.sequence + sequence
                
                    if neighbourConsensusEntry is not None:
                        neighbourConsensusEntry.sequence = ("-"*nrGaps) + neighbourConsensusEntry.sequence
                    
            # entry should not be merged or could be neither appended nor prepended 
            # add LCBs to alignment as it is
            else:
                mergedSplitLCBs.append(lcb)
                

        merged = Alignment(alignment.xmfaFile)
        for nr, genome in alignment.genomes.items():
            merged.addGenome(genome, nr)
        
        for lcb in mergedSplitLCBs:
            merged.addLCB(lcb)
                
        return merged
      

class Realigner:

    _sub_matrix = {  ('A', 'A'): 5,
                     ('C', 'C'): 5, ('C', 'A'): -4,
                     ('G', 'G'): 5, ('G', 'A'): -4, ('G', 'C'): -4,
                     ('T', 'T'): 5, ('T', 'A'): -4, ('T', 'C'): -4, ('T', 'G'): -4,
                     ('M', 'M'): -1, ('M', 'A'): 1, ('M', 'C'): 1, ('M', 'G'): -4, ('M', 'T'): -4,
                     ('R', 'R'): -1, ('R', 'A'): 1, ('R', 'C'): -4, ('R', 'G'): 1, ('R', 'T'): -4, ('R', 'M'): -2, 
                     ('W', 'W'): -1, ('W', 'A'): 1, ('W', 'C'): -4, ('W', 'G'): -4, ('W', 'T'): 1, ('W', 'M'): -2, ('W', 'R'): -2, 
                     ('S', 'S'): -1, ('S', 'A'): -4, ('S', 'C'): 1, ('S', 'G'): 1, ('S', 'T'): -4, ('S', 'M'): -2, ('S', 'R'): -2, ('S', 'W'): -4, 
                     ('Y', 'Y'): -1, ('Y', 'A'): -4, ('Y', 'C'): 1, ('Y', 'G'): -4, ('Y', 'T'): 1, ('Y', 'M'): -2, ('Y', 'R'): -4, ('Y', 'W'): -2, ('Y', 'S'): -2, 
                     ('K', 'K'): -1, ('K', 'A'): -4, ('K', 'C'): -4, ('K', 'G'): 1, ('K', 'T'): 1, ('K', 'M'): -4, ('K', 'R'): -2, ('K', 'W'): -2, ('K', 'S'): -2, ('K', 'Y'): -2, 
                     ('V', 'V'): -1, ('V', 'A'): -1, ('V', 'C'): -1, ('V', 'G'): -1, ('V', 'T'): -4, ('V', 'M'): -1, ('V', 'R'): -1, ('V', 'W'): -3, ('V', 'S'): -1, ('V', 'Y'): -3, ('V', 'K'): -3, 
                     ('H', 'H'): -1, ('H', 'A'): -1, ('H', 'C'): -1, ('H', 'G'): -4, ('H', 'T'): -1, ('H', 'M'): -1, ('H', 'R'): -3, ('H', 'W'): -1, ('H', 'S'): -3, ('H', 'Y'): -1, ('H', 'K'): -3, ('H', 'V'): -2, 
                     ('D', 'D'): -1, ('D', 'A'): -1, ('D', 'C'): -4, ('D', 'G'): -1, ('D', 'T'): -1, ('D', 'M'): -3, ('D', 'R'): -1, ('D', 'W'): -1, ('D', 'S'): -3, ('D', 'Y'): -3, ('D', 'K'): -1, ('D', 'V'): -2, ('D', 'H'): -2, 
                     ('B', 'B'): -1, ('B', 'A'): -4, ('B', 'C'): -1, ('B', 'G'): -1, ('B', 'T'): -1, ('B', 'M'): -3, ('B', 'R'): -3, ('B', 'W'): -3, ('B', 'S'): -1, ('B', 'Y'): -1, ('B', 'K'): -1, ('B', 'V'): -2, ('B', 'H'): -2, ('B', 'D'): -2, 
                     ('N', 'N'): -1, ('N', 'A'): -2, ('N', 'C'): -2, ('N', 'G'): -2, ('N', 'T'): -2, ('N', 'M'): -1, ('N', 'R'): -1, ('N', 'W'): -1, ('N', 'S'): -1, ('N', 'Y'): -1, ('N', 'K'): -1, ('N', 'V'): -1, ('N', 'H'): -1, ('N', 'D'): -1, ('N', 'B'): -1
                    }
                    
    ## local realignment around overlapping or consecutive gaps in two sequences
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
                oneFirstTwoSecond = self._getRealignRegions(entryOne.gaps, entryTwo.gaps)
                
                if len(oneFirstTwoSecond) > 0:
                    seqOne, seqTwo = self._realign(entryOne.sequence, entryTwo.sequence, oneFirstTwoSecond)
                    entryOne.sequence = seqOne
                    entryTwo.sequence = seqTwo
                
                # get regions to realign for updated entries with second entry first
                twoFirstOneSecond = self._getRealignRegions(entryTwo.gaps, entryOne.gaps)
                
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

                minSeqLength = minSeqLength*2
            
                seqStart = interval[minIndex][0] - index_offset - minSeqLength
                seqEnd = interval[minIndex][1] - index_offset + minSeqLength
                
                # do not go over boundaries of sequences!
                seqStart = max(seqStart, 0)
                minOrgSeqLength = min(len(seqOne), len(seqTwo)) - 1
                seqEnd = min(seqEnd, minOrgSeqLength)
                                
                alignments = pairwise2.align.globalds(seqOne[seqStart:seqEnd].replace("-", "").upper() , seqTwo[seqStart:seqEnd].replace("-", "").upper(), self._sub_matrix, -0.5, -0.1, one_alignment_only=True)
                if len(alignments) > 0:                                
                    maxscore = max( [x[2] for x in alignments] )
                    alignments = (lambda maxscore=maxscore: [item for item in alignments if item[2] == maxscore])()
                
                    minlength = min( [x[4] for x in alignments] )
                    alignments = (lambda minlength=minlength: [item for item in alignments if item[4] == minlength])()
                
                    seqOne = alignments[0][0].join([ seqOne[:seqStart], seqOne[seqEnd:] ])
                    seqTwo = alignments[0][1].join([ seqTwo[:seqStart], seqTwo[seqEnd:] ])
                
                    index_offset += ((seqEnd - seqStart) - minlength)

        return (seqOne, seqTwo)
