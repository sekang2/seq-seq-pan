import bisect

from seqseqpan.base import *
from seqseqpan.resolver import Resolver


class Splitter:
    def __init__(self, alignment, chromosome_desc):
        self.alignment = alignment
        self.resolver = Resolver()
        for nr, genome in alignment.genomes.items():
            genome.add_chromosomes(chromosome_desc[nr])

    def split_alignment(self):

        splitted_lcbs = []

        for lcb in self.alignment.lcbs:
            splitted_lcbs.extend(self.split_by_chromosomes(lcb))

        self.alignment.lcbs = splitted_lcbs

        return self.alignment

    def split_by_chromosomes(self, lcb):

        split_coords = []

        for entry in lcb.entries:
            chrom_starts = self.get_chromosomes_for_entry(entry)
            starts = chrom_starts

            if entry.strand == "-":
                starts = [(entry.end - chrom_start) + 2 for chrom_start in starts]
            else:
                starts = [(chrom_start - entry.start) + 1 for chrom_start in starts]

            starts[0] = 1

            if len(starts) > 1:

                starts_with_gaps = [entry.get_position_within_entry_with_gaps(coords_local) for coords_local in starts]
                starts_with_gaps[0] = 1
                split_coords.extend( starts_with_gaps )

        split_coords = sorted(list(set(split_coords)))

        if len(split_coords) > 1:

            split_coords = [c - 1 for c in split_coords]
            split_coords = split_coords + [lcb.length]

            # create new lcbs with split coords and return them
            lcbs = [LCB() for _ in range(len(split_coords) - 1)]
            entries = lcb.entries

            for c_idx in range(1, len(split_coords)):
                for e_idx in range(len(entries)):

                    entry = lcb.entries[e_idx]

                    start = split_coords[c_idx - 1]
                    end = split_coords[c_idx]

                    new_entry = self.resolver.get_split_entry(entry, start, end)
                    if new_entry is not None:
                        lcbs[c_idx - 1].add_entries(new_entry)
            return lcbs
        else:
            return [lcb]

    def get_chromosomes_for_entry(self, entry):
        genome = self.alignment.genomes[entry.genome_nr]
        chrom_starts = sorted(genome.chromosomes.keys())
        idx_1 = bisect.bisect_right(chrom_starts, entry.start)
        idx_1 -= 1
        idx_2 = bisect.bisect_right(chrom_starts, entry.end)

        return chrom_starts[idx_1:idx_2]
