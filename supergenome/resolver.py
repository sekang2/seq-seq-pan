import re
import bisect

from supergenome.exception import *
from supergenome.base import *

class Resolver:
    
    def resolveMultiAlignment(self, alignment, consensus):
        
        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()

        if alignment.genomes[1].filepath == consensus.fastaFile:
            consensusGenomeNr = 1
        elif alignment.genomes[2].filepath == consensus.fastaFile:
            consensusGenomeNr = 2
        else:
            raise ConsensusGenomeNumberError()
        
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
                splitLcbs = self._splitLCB(lcb, consensusGenomeNr, newGenomeNr, consensus) # split LCB                    
                for slcb in splitLcbs:
                    resolved.addLCB(slcb)
            else:
                resolved.addLCB(lcb)
        
        return(resolved)
        
            
    def reconstructAlignment(self, resolvedaln, consensus, orgAlignment):
        if len(resolvedaln.genomes) > 2:
            raise ConsensusXMFAInputError()

        if resolvedaln.genomes[1].filepath == consensus.fastaFile:
            consensusGenomeNr = 1
        elif resolvedaln.genomes[2].filepath == consensus.fastaFile:
            consensusGenomeNr = 2
        else:
            raise ConsensusGenomeNumberError()
        
        newGenomeNr = (1 if consensusGenomeNr == 2 else 2)
        
        # add genomes from consensus and new genome
        recalculated = Alignment(orgAlignment.xmfaFile)
        for nr, genome in orgAlignment.genomes.items():
            recalculated.addGenome(genome, nr)
        nrGenomes = len(orgAlignment.genomes)+1
        recalculated.addGenome(resolvedaln.genomes[newGenomeNr], nrGenomes)
        try:
            sortedOrgLCBs = orgAlignment.getSortedLCBs(consensus.order)    
        except ParameterError:
            raise ConsensusFastaFormatError()
        else:
        
            # for each block reconstruct original sequence from consensus entry and/or add entry from new genome
            for lcb in resolvedaln.LCBs:
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
                    except IndexError as e:
                        print(lcb.number)
                        raise e
                
                recalculated.addLCB(recalculatedLCB)
            
            return recalculated
        
        
    def _splitLCB(self, lcb, consensusGenomeNr, newGenomeNr, consensus):
        splitLCBs = []
        delimiterPositions = []
             
        consensusEntry = lcb.getEntry(consensusGenomeNr)
        
        ## check whether consensusEntry overlaps delimiter position      
        
        # check if entry contains delimiter sequence less than normal length at beginning
        intervals = [idx - 1000 for idx in sorted(consensus.blockStartIndices)[1:] if 0 <= (idx - consensusEntry.start) <= len(BLOCK_DELIMITER)]
        
        # check for all overlaps
        intervals.extend([idx - 1000 for idx in sorted(consensus.blockStartIndices)[1:] if consensusEntry.start < (idx - 1000) <  consensusEntry.end])
        
        for startInterval in intervals:
            startWithinBlock = startInterval - consensusEntry.start + 1
            offset = (startWithinBlock * -1 if startWithinBlock < 0 else 0)
            # in case of delimiter sequencing at beginning of entry, index will be negative but should be 0
            startWithinBlock = (consensusEntry.getPositionWithinEntryWithGaps(startWithinBlock) if startWithinBlock > 0 else 0)
                
            # search for delimiter sequence starting at interval position: N with gaps inbetween with maximum length of 1000
            m = re.search("(N-*){0,"+str(len(BLOCK_DELIMITER)-offset-1)+"}N", consensusEntry.sequence[startWithinBlock:])
            if m is None:
                raise ConsensusCorruptError(startWithinBlock+consensusEntry.start)
            
            # store positions to split
            delimiterPositions.append(startWithinBlock + m.start(0))
            delimiterPositions.append(startWithinBlock + m.end(0))
        
        
        if len(delimiterPositions) == 0: # no entries with delimiter sequence in LCB
            return [lcb]
        
        
        indices = delimiterPositions
        
        # store whether even or uneven elements are delimiters
        evenIsDelimiter = True
        if delimiterPositions[0] > 0:
            indices = [0] + indices
            evenIsDelimiter = False
        
        if delimiterPositions[-1] < lcb.length:
            indices = indices + [lcb.length]
        
        for i in range(len(indices)-1):
        
            newEntry = lcb.getEntry(newGenomeNr)
            if newEntry is None:
                entries = []
            else:
                entries = [self._getNewEntry(newEntry, indices[i], indices[i+1])]
            
            # check if element is non-delimiter element
            if (i % 2 == 0 and not evenIsDelimiter) or ( i % 2 != 0 and evenIsDelimiter):
                entries.append(self._getNewEntry(consensusEntry, indices[i], indices[i+1]))
            
            entries = list(filter(None, entries)) # add all entries that are not None)
            
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
        if seq == "" or re.search("[^-]", seq) is None:
            return None
        splitlen = splitend - splitstart
        
        sumgaps = 0
        
        subgaps = entry.getSubGapList(0, splitstart - 1)
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
            eSubgaps = e.getSubGapList(0,startWithinBlock - 1)
            if len(eSubgaps) > 0:
                eSumgaps = sum(end-start for start, end in eSubgaps.items())
            
            start = start - eSumgaps
            end = end - eSumgaps
            
            newEntry = SequenceEntry(e.genomeNr, start, end, e.strand, sequence)
            
            newEntry.end = newEntry.end - sum( end-start for start, end in newEntry.gaps.items() )
            
            # add all gaps in consensus sequence to original sequence for correct alignment
            for gstart, gend in sorted(consensusEntry.gaps.items()):
                sequence = self._insertGap(sequence, gstart, gend-gstart)
            
            
            if sequence != "" and re.search("[^-]", sequence) is not None:
                newEntry.sequence = sequence
                orgEntries.append(newEntry)
        
        return orgEntries
    
    
    def _insertGap(self, sequence, position, length):
        seq = ("-"*length).join([sequence[:position], sequence[position:]])
        return seq
