import bisect

from seqseqpan.base import *


class Splitter:
    def __init__(self, alignment, chromosome_desc):
        self.alignment = alignment
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
            starts[0] = entry.start
            if len(starts) > 1:
                split_coords.extend(
                    [entry.get_position_within_entry_with_gaps((chrom_start - entry.start) + 1) for chrom_start in chrom_starts])

        split_coords = sorted(list(set(split_coords)))

        if len(split_coords) > 1:
            print("Splitting LCB by chromosomes: " + str(lcb.number) + " ")
            print(split_coords)
            print("\n")
            split_coords = [c - 1 for c in split_coords]
            split_coords = split_coords + [lcb.length]

            # create new lcbs with split coords and return them
            lcbs = [LCB() for _ in range(len(split_coords) - 1)]
            entries = lcb.entries
            cur_starts = [entry.start for entry in entries]

            for c_idx in range(1, len(split_coords)):
                for e_idx in range(len(entries)):

                    entry = lcb.entries[e_idx]
                    start = cur_starts[e_idx]

                    seq = entry.sequence[split_coords[c_idx - 1]:split_coords[c_idx]]

                    non_gaps = len(seq) - seq.count("-")

                    if non_gaps > 0:
                        cur_starts[e_idx] = start + non_gaps
                        end = cur_starts[e_idx] - 1
                        new_entry = SequenceEntry(entry.genome_nr, start, end, entry.strand, seq)
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
