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
            
            # check if new entry is small and only created by splitting or aligning of new genome (consensus is None)
            if consensusEntry is None and newEntry is not None and len(newEntry.sequence) <= 10:
                
                nrGaps = len(newEntry.sequence)
                #sequence = newEntry.sequence
                
                # try to append to previous entry
                if i > 0:
                    lastLCB = mergedSplitLCBs[-1]
                    lastNewEntry = lastLCB.getEntry(newGenomeNr)
                    if lastNewEntry is not None:
                        tryNextEntry = False
                        
                        sequence = lastNewEntry.sequence
                        
                        lastNewEntry.sequence = sequence + newEntry.sequence
                        lastNewEntry.end = newEntry.end
                    
                        lastConsensusEntry = lastLCB.getEntry(consensusGenomeNr)
                        if lastConsensusEntry is not None:
                            lastConsensusEntry.sequence += ("-"*nrGaps)
                    
                    lastLCB.length += nrGaps
                    
                # previous entry did not work, try to prepend to next entry
                if tryNextEntry and len(lcbs) > i:
                    nextLCB = lcbs[i+1]
                    nextNewEntry = nextLCB.getEntry(newGenomeNr)
                    if nextNewEntry is not None:
                        tryNextEntry = False
                        
                        sequence = nextNewEntry.sequence
                        
                        nextNewEntry.sequence = newEntry.sequence + sequence
                        nextNewEntry.start = newEntry.start
                    
                        nextConsensusEntry = nextLCB.getEntry(consensusGenomeNr)
                        if nextConsensusEntry is not None:
                            nextConsensusEntry.sequence = ("-"*nrGaps) + nextConsensusEntry.sequence
                        
                        nextLCB.length += nrGaps
                        
            # entry should not be merged or could be neither appended nor prepended 
            # add LCBs to alignment as it is
            if tryNextEntry:
                mergedSplitLCBs.append(lcb)

        merged = Alignment(alignment.xmfaFile)
        for nr, genome in alignment.genomes.items():
            merged.addGenome(genome, nr)
        
        for lcb in mergedSplitLCBs:
            merged.addLCB(lcb)
                
        return merged
      

class Realigner:
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
            
                seqStart = interval[minIndex][0] - index_offset - minSeqLength
                seqEnd = interval[minIndex][1] - index_offset + minSeqLength
                
                # do not go over boundaries of sequences!
                seqStart = max(seqStart, 0)
                minOrgSeqLength = min(len(seqOne), len(seqTwo)) - 1
                seqEnd = min(seqEnd, minOrgSeqLength)
                
                alignments = pairwise2.align.globalxx(seqOne[seqStart:seqEnd].replace("-", "") , seqTwo[seqStart:seqEnd].replace("-", ""))
                                
                maxscore = max( [x[2] for x in alignments] )
                alignments = (lambda maxscore=maxscore: [item for item in alignments if item[2] == maxscore])()
                
                minlength = min( [x[4] for x in alignments] )
                alignments = (lambda minlength=minlength: [item for item in alignments if item[4] == minlength])()
                
                seqOne = alignments[0][0].join([ seqOne[:seqStart], seqOne[seqEnd:] ])
                seqTwo = alignments[0][1].join([ seqTwo[:seqStart], seqTwo[seqEnd:] ])
                
                index_offset += ((seqEnd - seqStart) - minlength)

        return (seqOne, seqTwo)