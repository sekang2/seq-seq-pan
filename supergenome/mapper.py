import bisect
from collections import defaultdict
import pdb
    
class Mapper:
    
    def mapCoordinates(self, alignment, consensus, source, dests, coordinates):
        
        print("len(coordinates):" + str(len(coordinates)))
        
        if type(dests) is not list:
            dests = [dests]

        dests = sorted(dests)
        coordinates = sorted(coordinates)
        coord_dict = defaultdict(dict)
        
        addC = ("c" in dests)
        
        # store and do not reorder!
        lcbs = alignment.LCBs
        
        # remove consensus-dest from list (calculation is different)
        if addC:
            dests = sorted(dests)[:-1]
        
        if source == "c":
            i = 0
            for coord in coordinates:
                if i%1000 == 0:
                    print(i)
                idx = bisect.bisect_left(consensus.blockStartIndices, coord) 
                idx -= 1
                lcb = lcbs[idx]
                
                # calculate position within current consensus block - there are no gaps in consensus sequence!
                posWithinBlock = coord - consensus.blockStartIndices[idx] - 1
                
                if addC:
                    coord_dict[coord]["c"] = coord
                
                coord_dict[coord].update(self._getCoordsForEntries(lcb.entries, dests, posWithinBlock))
                i = i+1            
        else:
            
            sourceBlocks = {}
            # store ends of blocks for finding lcb for coords
            
            for i in range(len(lcbs)):
                e = lcbs[i].getEntry(int(source))
                if e is not None:
                    sourceBlocks[e.end] = {"lcb":i, "entry":e }
            
            for coord in coordinates:
                sourceEnds = sorted(sourceBlocks.keys())
                idxInEnds = bisect.bisect_left(sourceEnds, coord)
                lcbIdx = sourceBlocks[sourceEnds[idxInEnds]]["lcb"]
                
                if lcbIdx > len(lcbs):
                    raise CoordinateOutOfBoundsError(coord, source)
    
                sourceEntry = sourceBlocks[sourceEnds[idxInEnds]]["entry"]
                lcb = lcbs[lcbIdx]
                
                if sourceEntry.strand == "+":
                    posWithinBlockWithoutGaps = coord - sourceEntry.start
                else:
                    posWithinBlockWithoutGaps = sourceEntry.end - coord
            
                # add consensus coordinates to dict
                if addC:
                    consLength = sum([lcb.length for lcb in lcbs[0:lcbIdx]])
                    coord_dict[coord]["c"] = consLength + posWithinBlockWithoutGaps + 1
                    
                # check if dests other than consensus are needed
                if len(dests) > 0:
                
                    posWithinBlock = sourceEntry.getPositionWithinEntryWithGaps(posWithinBlockWithoutGaps)
                    
                    coord_dict[coord].update(self._getCoordsForEntries(lcb.entries, dests, posWithinBlock))
            

        return coord_dict
        
        
    def _getCoordsForEntries(self, entries, dests, posWithinBlock):
        coord_dict = {}
        
        for e in entries: # faster than looping through entries everytime to get entry of genome x
            if str(e.genomeNr) in dests: 
            
                isGap = False
            
                if e.strand == "+":
                    eSubgaps = e.getSubGapList(0,posWithinBlock)
                    eSumgaps = 0
                    
                    if len(eSubgaps) > 0:
                        eSumgaps = sum(end-start for start, end in eSubgaps.items())
                        lastGap = sorted(eSubgaps.items())[-1]
                        isGap = (lastGap[0] <= posWithinBlock and posWithinBlock < lastGap[1])
                    
                    if not isGap:
                        coord_dict[str(e.genomeNr)] = e.start + (posWithinBlock - eSumgaps)
                else:
                    eSubgaps = e.getSubGapList(e.end - posWithinBlock, e.end)
                    eSumgaps = 0
                    if len(eSubgaps) > 0:
                        eSumgaps = sum(end-start for start, end in eSubgaps.items())
                        lastGap = sorted(eSubgaps.items())[-1]
                        isGap = (lastGap[0] <= posWithinBlock and posWithinBlock < lastGap[1])
                    
                    if not isGap:
                        coord_dict[str(e.genomeNr)] = (e.end - (posWithinBlock - eSumgaps)) * -1
        
        return coord_dict

        