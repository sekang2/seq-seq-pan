import bisect

from supergenome.base import *

class Splitter:

    def __init__(self, alignment):
        self.alignment = alignment
        for nr, genome in alignment.genomes.items():
            genome.readChromosomes()

            
    def splitByChromosomes(self, lcb):
        
        split_coords = []
        
        for entry in lcb.entries:
            chromosomeStarts = self.getChromosomesForEntry(entry)
            starts = chromosomeStarts
            starts[0] = entry.start
            split_coords.extend([entry.getPositionWithinEntryWithGaps((chrstart - entry.start) + 1) for chrstart in chromosomeStarts])
        
        split_coords = sorted(list(set(split_coords)))
        
        
        if len(split_coords) > 1:
            print("Splitting LCB by chromosomes: " + str(lcb.number) + " ")
            print(split_coords)
            print("\n")
            split_coords = split_coords + [lcb.length]
            split_coords = [c - 1 for c in split_coords]
            
            # create new lcbs with split coords and return them
            lcbs = [LCB() for _ in range(len(split_coords)-1)]
            entries = lcb.entries
            cur_starts = [entry.start for entry in entries]

            for c_idx in range(1,len(split_coords)):
                for e_idx in range(len(entries)):
                    
                    entry = lcb.entries[e_idx]
                    start = cur_starts[e_idx]
                    
                    seq = entry.sequence[split_coords[c_idx - 1]:split_coords[c_idx]]
                    
                    non_gaps = len(seq) - seq.count("-")
                    
                    if non_gaps > 0:
                        cur_starts[e_idx] = start + non_gaps
                        end = cur_starts[e_idx] - 1
                        newEntry = SequenceEntry(entry.genomeNr, start, end, entry.strand, seq)
                        lcbs[c_idx - 1].addEntries(newEntry)
            return lcbs
        else:
            return [lcb] 

            
    def getChromosomesForEntry(self, entry):
        genome = self.alignment.genomes[entry.genomeNr]
        chrstarts = sorted(genome.chromosomes.keys())
        idx_1 = bisect.bisect_right(chrstarts, entry.start)
        idx_1 -= 1
        idx_2 = bisect.bisect_right(chrstarts, entry.end)
        
        return chrstarts[idx_1:idx_2]
