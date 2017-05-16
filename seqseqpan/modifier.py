import collections
import itertools
from operator import itemgetter
import sys

from Bio import pairwise2

from skbio.alignment import local_pairwise_align_ssw, global_pairwise_align
from skbio import DNA

from seqseqpan.exception import ConsensusXMFAInputError
from seqseqpan.base import *
import pdb

class Separator:
    def separate_lcbs(self, alignment, length):
        separated = Alignment(alignment.xmfa_file)
        for nr, genome in alignment.genomes.items():
            separated.add_genome(genome, nr)

        for lcb in alignment.lcbs:
            if lcb.length <= length:
                for entry in lcb.entries:
                    seq = entry.sequence.replace("-", "")
                    new_entry = SequenceEntry(entry.genome_nr, entry.start, entry.end, entry.strand, seq)
                    separated.add_lcb_entries(new_entry)
            else:
                separated.add_lcb(lcb)

        return separated


class Merger:
    def merge_lcbs(self, alignment, consensus_genome_nr, new_genome_nr, block_length):

        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()

        lcbs = alignment.get_sorted_lcbs(new_genome_nr)

        # do not create small (less bp than blocklength) LCBs by splitting, but append/prepend sequence         
        merged_split_lcbs = []

        for i in range(len(lcbs)):

            lcb = lcbs[i]
            new_entry = lcb.get_entry(new_genome_nr)
            consensus_entry = lcb.get_entry(consensus_genome_nr)
            prepend = False
            append = False
            to_reverse_complement = False

            neighbour_new_entry = None

            # check if new entry is small and only created by splitting or aligning of new genome (consensus is None)
            if consensus_entry is None and new_entry is not None and len(new_entry.sequence) <= block_length:

                nr_gaps = len(new_entry.sequence)

                # try to append to previous entry
                if i > 0:
                    neighbour_lcb = merged_split_lcbs[-1]
                    neighbour_new_entry = neighbour_lcb.get_entry(new_genome_nr)
                    if neighbour_new_entry is not None and (new_entry.start - neighbour_new_entry.end) == 1:

                        if neighbour_new_entry.strand == "+":
                            append = True
                        else:
                            prepend = True

                        if new_entry.strand != neighbour_new_entry.strand:
                            to_reverse_complement = True

                # previous entry did not work, try to prepend to next entry
                if not prepend and not append and len(lcbs) > i:
                    neighbour_lcb = lcbs[i + 1]
                    neighbour_new_entry = neighbour_lcb.get_entry(new_genome_nr)
                    if neighbour_new_entry is not None and (neighbour_new_entry.start - new_entry.end) == 1:

                        if neighbour_new_entry.strand == "+":
                            prepend = True
                        else:
                            append = True

                        if new_entry.strand != neighbour_new_entry.strand:
                            to_reverse_complement = True

            if append or prepend:

                if to_reverse_complement:
                    new_entry.reverse_complement()

                sequence = neighbour_new_entry.sequence
                neighbour_consensus_entry = neighbour_lcb.get_entry(consensus_genome_nr)
                neighbour_lcb.length += nr_gaps

                neighbour_new_entry.end = max(new_entry.end, neighbour_new_entry.end)
                neighbour_new_entry.start = min(new_entry.start, neighbour_new_entry.start)

                if append:

                    neighbour_new_entry.sequence = sequence + new_entry.sequence

                    if neighbour_consensus_entry is not None:
                        neighbour_consensus_entry.sequence += ("-" * nr_gaps)

                elif prepend:

                    neighbour_new_entry.sequence = new_entry.sequence + sequence

                    if neighbour_consensus_entry is not None:
                        neighbour_consensus_entry.sequence = ("-" * nr_gaps) + neighbour_consensus_entry.sequence

            # entry should not be merged or could be neither appended nor prepended 
            # add LCBs to alignment as it is
            else:
                merged_split_lcbs.append(lcb)

        merged = Alignment(alignment.xmfa_file)
        for nr, genome in alignment.genomes.items():
            merged.add_genome(genome, nr)

        for lcb in merged_split_lcbs:
            merged.add_lcb(lcb)

        return merged


class Realigner:
    _sub_matrix = {('A', 'A'): 5,
                   ('C', 'C'): 5, ('C', 'A'): -4,
                   ('G', 'G'): 5, ('G', 'A'): -4, ('G', 'C'): -4,
                   ('T', 'T'): 5, ('T', 'A'): -4, ('T', 'C'): -4, ('T', 'G'): -4,
                   ('M', 'M'): -1, ('M', 'A'): 1, ('M', 'C'): 1, ('M', 'G'): -4, ('M', 'T'): -4,
                   ('R', 'R'): -1, ('R', 'A'): 1, ('R', 'C'): -4, ('R', 'G'): 1, ('R', 'T'): -4, ('R', 'M'): -2,
                   ('W', 'W'): -1, ('W', 'A'): 1, ('W', 'C'): -4, ('W', 'G'): -4, ('W', 'T'): 1, ('W', 'M'): -2,
                   ('W', 'R'): -2,
                   ('S', 'S'): -1, ('S', 'A'): -4, ('S', 'C'): 1, ('S', 'G'): 1, ('S', 'T'): -4, ('S', 'M'): -2,
                   ('S', 'R'): -2, ('S', 'W'): -4,
                   ('Y', 'Y'): -1, ('Y', 'A'): -4, ('Y', 'C'): 1, ('Y', 'G'): -4, ('Y', 'T'): 1, ('Y', 'M'): -2,
                   ('Y', 'R'): -4, ('Y', 'W'): -2, ('Y', 'S'): -2,
                   ('K', 'K'): -1, ('K', 'A'): -4, ('K', 'C'): -4, ('K', 'G'): 1, ('K', 'T'): 1, ('K', 'M'): -4,
                   ('K', 'R'): -2, ('K', 'W'): -2, ('K', 'S'): -2, ('K', 'Y'): -2,
                   ('V', 'V'): -1, ('V', 'A'): -1, ('V', 'C'): -1, ('V', 'G'): -1, ('V', 'T'): -4, ('V', 'M'): -1,
                   ('V', 'R'): -1, ('V', 'W'): -3, ('V', 'S'): -1, ('V', 'Y'): -3, ('V', 'K'): -3,
                   ('H', 'H'): -1, ('H', 'A'): -1, ('H', 'C'): -1, ('H', 'G'): -4, ('H', 'T'): -1, ('H', 'M'): -1,
                   ('H', 'R'): -3, ('H', 'W'): -1, ('H', 'S'): -3, ('H', 'Y'): -1, ('H', 'K'): -3, ('H', 'V'): -2,
                   ('D', 'D'): -1, ('D', 'A'): -1, ('D', 'C'): -4, ('D', 'G'): -1, ('D', 'T'): -1, ('D', 'M'): -3,
                   ('D', 'R'): -1, ('D', 'W'): -1, ('D', 'S'): -3, ('D', 'Y'): -3, ('D', 'K'): -1, ('D', 'V'): -2,
                   ('D', 'H'): -2,
                   ('B', 'B'): -1, ('B', 'A'): -4, ('B', 'C'): -1, ('B', 'G'): -1, ('B', 'T'): -1, ('B', 'M'): -3,
                   ('B', 'R'): -3, ('B', 'W'): -3, ('B', 'S'): -1, ('B', 'Y'): -1, ('B', 'K'): -1, ('B', 'V'): -2,
                   ('B', 'H'): -2, ('B', 'D'): -2,
                   ('N', 'N'): -1, ('N', 'A'): -2, ('N', 'C'): -2, ('N', 'G'): -2, ('N', 'T'): -2, ('N', 'M'): -1,
                   ('N', 'R'): -1, ('N', 'W'): -1, ('N', 'S'): -1, ('N', 'Y'): -1, ('N', 'K'): -1, ('N', 'V'): -1,
                   ('N', 'H'): -1, ('N', 'D'): -1, ('N', 'B'): -1
                   }

    _sub_matrix_dict = {
        'A': {'A': 5, 'C': -4, 'G': -4, 'T': -4, 'M': 1, 'R': 1, 'W': 1, 'S': -4, 'Y': -4, 'K': -4, 'V': -1, 'H': -1,
              'D': -1, 'B': -4, 'N': -1},
        'C': {'A': -4, 'C': 5, 'G': -4, 'T': -4, 'M': 1, 'R': -4, 'W': -4, 'S': 1, 'Y': 1, 'K': -4, 'V': -1, 'H': -1,
              'D': -4, 'B': -1, 'N': -1},
        'G': {'A': -4, 'C': -4, 'G': 5, 'T': -4, 'M': -4, 'R': 1, 'W': -4, 'S': 1, 'Y': -4, 'K': 1, 'V': -1, 'H': -4,
              'D': -1, 'B': -1, 'N': -1},
        'T': {'A': -4, 'C': -4, 'G': -4, 'T': 5, 'M': -4, 'R': -4, 'W': 1, 'S': -4, 'Y': 1, 'K': 1, 'V': -4, 'H': -1,
              'D': -1, 'B': -1, 'N': -1},
        'M': {'A': 1, 'C': 1, 'G': -4, 'T': -4, 'M': -1, 'R': -2, 'W': -2, 'S': -2, 'Y': -2, 'K': -4, 'V': -1, 'H': -1,
              'D': -3, 'B': -3, 'N': -1},
        'R': {'A': 1, 'C': -4, 'G': 1, 'T': -4, 'M': -2, 'R': -1, 'W': -2, 'S': -2, 'Y': -4, 'K': -2, 'V': -1, 'H': -3,
              'D': -1, 'B': -3, 'N': -1},
        'W': {'A': 1, 'C': -4, 'G': -2, 'T': -4, 'M': -2, 'R': -2, 'W': -1, 'S': -4, 'Y': -2, 'K': -2, 'V': -3, 'H': -1,
              'D': -1, 'B': -3, 'N': -1},
        'S': {'A': -4, 'C': 1, 'G': 1, 'T': -4, 'M': -2, 'R': -2, 'W': -4, 'S': -1, 'Y': -2, 'K': -2, 'V': -1, 'H': -3,
              'D': -3, 'B': -1, 'N': -1},
        'Y': {'A': -4, 'C': 1, 'G': -4, 'T': 1, 'M': -2, 'R': -4, 'W': -2, 'S': -2, 'Y': -1, 'K': -2, 'V': -3, 'H': -1,
              'D': -3, 'B': -1, 'N': -1},
        'K': {'A': -4, 'C': -4, 'G': 1, 'T': 1, 'M': -4, 'R': -2, 'W': -2, 'S': -2, 'Y': -2, 'K': -1, 'V': -3, 'H': -3,
              'D': -1, 'B': -1, 'N': -1},
        'V': {'A': -1, 'C': -1, 'G': -1, 'T': -4, 'M': -1, 'R': -1, 'W': -3, 'S': -1, 'Y': -3, 'K': -3, 'V': -1,
              'H': -2, 'D': -2, 'B': -2, 'N': -1},
        'H': {'A': -1, 'C': -1, 'G': -4, 'T': -1, 'M': -1, 'R': -3, 'W': -1, 'S': -3, 'Y': -1, 'K': -3, 'V': -2,
              'H': -1, 'D': -2, 'B': -2, 'N': -1},
        'D': {'A': -1, 'C': -4, 'G': -1, 'T': -1, 'M': -3, 'R': -1, 'W': -1, 'S': -3, 'Y': -3, 'K': -1, 'V': -1,
              'H': -2, 'D': -1, 'B': -2, 'N': -1},
        'B': {'A': -4, 'C': -1, 'G': -1, 'T': -1, 'M': -3, 'R': -3, 'W': -3, 'S': -1, 'Y': -1, 'K': -1, 'V': -2,
              'H': -2, 'D': -2, 'B': -1, 'N': -1},
        'N': {'A': -2, 'C': -2, 'G': -2, 'T': -2, 'M': -1, 'R': -1, 'W': -1, 'S': -1, 'Y': -1, 'K': -1, 'V': -1,
              'H': -1, 'D': -1, 'B': -1, 'N': -1}}


    # local realignment around overlapping or consecutive gaps in two sequences
    def realign(self, alignment):
        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()

        realigned = Alignment(alignment.xmfa_file)
        for nr, genome in alignment.genomes.items():
            realigned.add_genome(genome, nr)

        # go through lcbs, skip one-entry ones
        for lcb in alignment.get_sorted_lcbs(0):

            if len(lcb.entries) == 1:
                realigned.add_lcb(lcb)
            else:
                entry_one = lcb.entries[0]
                entry_two = lcb.entries[1]

                # get regions to realign
                one_first_two_second = self._get_realign_regions(entry_one.gaps, entry_two.gaps)

                if len(one_first_two_second) > 0:
                    seq_one, seq_two = self._realign(entry_one.sequence, entry_two.sequence, one_first_two_second)
                    entry_one.sequence = seq_one
                    entry_two.sequence = seq_two

                # get regions to realign for updated entries with second entry first
                two_first_one_second = self._get_realign_regions(entry_two.gaps, entry_one.gaps)

                if len(two_first_one_second) > 0:
                    seq_two, seq_one = self._realign(entry_two.sequence, entry_one.sequence, two_first_one_second)
                    entry_one.sequence = seq_one
                    entry_two.sequence = seq_two

                new_lcb = LCB()
                new_lcb.add_entries([entry_one, entry_two])

                realigned.add_lcb(new_lcb)

        return realigned

    def _get_realign_regions(self, gaps_for_start, gaps_for_end):
        ends_dict = {end: start for start, end in gaps_for_end.items()}
        location_list = [start for start, end in gaps_for_start.items()] + list(ends_dict.keys())

        region_starts = [item for item, count in collections.Counter(location_list).items() if count > 1]

        regions = []
        for start in region_starts:
            regions.append([(start, gaps_for_start[start]), (ends_dict[start], start)])

        return regions

    def _realign(self, seq_one, seq_two, realign_regions):
        #pdb.set_trace()
        realign_regions = sorted(realign_regions)

        # offset for coordinates if sequences are changed due to alignment
        index_offset = 0

        for interval in realign_regions:
            #print(interval)

            #find start index and length of interval
            min_index, min_seq_length = min(enumerate([interval[0][1] - interval[0][0], interval[1][1] - interval[1][0]]),
                                            key=lambda p: p[1])

            #if min_seq_length < 16000:
            #print("\n")

            #print(min_seq_length)
            # double for more information
            min_seq_length *= 2

            # get sequence left and right of consecutive gap
            seq_start = interval[min_index][0] - index_offset - min_seq_length
            seq_end = interval[min_index][1] - index_offset + min_seq_length


            # do not go over boundaries of sequences!
            seq_start = max(seq_start, 0)
            min_orig_seq_length = min(len(seq_one), len(seq_two)) - 1
            seq_end = min(seq_end, min_orig_seq_length)

            #alignments = pairwise2.align.globalds(seq_one[seq_start:seq_end].replace("-", "").upper(),
            #                                      seq_two[seq_start:seq_end].replace("-", "").upper(),
            #                                      self._sub_matrix, -0.5, -0.1, one_alignment_only=True)

            #print(seq_one[seq_start:seq_end])
            #print(seq_two[seq_start:seq_end])

            seq_one_nogap = seq_one[seq_start:seq_end].replace("-", "")
            seq_two_nogap = seq_two[seq_start:seq_end].replace("-", "")

            #swap = False
            #if len(seq_two_nogap) > len(seq_one_nogap):
            #    swap = True
            #    seq_one_nogap, seq_two_nogap = seq_two_nogap, seq_one_nogap

            #print(seq_one_nogap)
            #print(seq_two_nogap)

            if not (seq_one_nogap == '' or seq_two_nogap == ''):

                try:
                    # fast local alignment
                    alignment, score, positions = local_pairwise_align_ssw(
                        DNA(seq_one_nogap.upper()),
                        DNA(seq_two_nogap.upper())#,
                        #substitution_matrix=self._sub_matrix_dict
                    )
                except IndexError as e:
                    pass
                else:
                    #locally aligned sequences
                    #if swap:
                    #aln_seq1 = str(alignment.iloc(axis='sequence')[1])
                    #aln_seq2 = str(alignment.iloc(axis='sequence')[0])

                    #aln_start1, aln_end1 = positions[1]
                    #aln_start2, aln_end2 = positions[0]

                    #seq_one_nogap, seq_two_nogap = seq_two_nogap, seq_one_nogap
                    #else:
                    aln_seq1 = str(alignment.iloc(axis='sequence')[0])
                    aln_seq2 = str(alignment.iloc(axis='sequence')[1])

                    aln_start1, aln_end1 = positions[0]
                    aln_start2, aln_end2 = positions[1]



                   # print(aln_seq1)
                   # print(aln_seq2)

                    # remove overlapping gaps in alignment (why are they even there?!?)
                    rm_indices = [i for i in range(len(aln_seq1)) if aln_seq1[i] == "-" and aln_seq2[i] == "-"]
                    aln_seq1 = "".join([char for idx, char in enumerate(aln_seq1) if idx not in rm_indices])
                    aln_seq2 = "".join([char for idx, char in enumerate(aln_seq2) if idx not in rm_indices])

                   # print(aln_seq1)
                   # print(aln_seq2)

                    # add sequences before and after local match to alignment

                    post_aln_len1 = len(seq_one_nogap) - aln_end1 - 1
                    post_aln_len2 = len(seq_two_nogap) - aln_end2 - 1

                    if aln_start2 > 0:
                        aln_seq1 = aln_start2*"-" + aln_seq1
                        aln_seq2 = seq_two_nogap[0:(aln_start2)] + aln_seq2

                    if aln_start1 > 0:
                        aln_seq1 = seq_one_nogap[0:(aln_start1)] + aln_seq1
                        aln_seq2 = aln_start1*"-" + aln_seq2

                    if post_aln_len1 > 0:
                        aln_seq1 += seq_one_nogap[(aln_end1+1):]
                        aln_seq2 += post_aln_len1*"-"

                    if post_aln_len2 > 0:
                        aln_seq1 += post_aln_len2*"-"
                        aln_seq2 += seq_two_nogap[(aln_end2+1):]


                   # print(aln_seq1)
                   # print(aln_seq2)

                   # print("\n\n")

                    # update entry sequences with better alignment
                    aln_length = len(aln_seq1)

                    seq_one = aln_seq1.join([seq_one[:seq_start], seq_one[seq_end:]])
                    seq_two = aln_seq2.join([seq_two[:seq_start], seq_two[seq_end:]])

                    index_offset += ((seq_end - seq_start) - aln_length)

        return seq_one, seq_two


class Remover:
    def remove(self, alignment, rm_genome):
        if len(alignment.genomes) >= rm_genome > -1:
            #print ("Number of LCBs: " + str(len(alignment.lcbs)))

            for lcb in alignment.lcbs:
                #sys.stdout.write(".")
                #sys.stdout.flush()

                entries = [entry for entry in lcb.entries if entry.genome_nr != rm_genome]

                # did LCB include entry of genome to remove?
                if len(entries) < len(lcb.entries):
                    # are there any entries left in LCB?
                    if len(entries) > 1:
                        # if more than one entry left search for gaps that are present in all remaining entries

                        rm_gaps = set(itertools.chain.from_iterable([list(range(k, entries[0].gaps[k])) for k in entries[0].gaps]))

                        for entry in entries[1:]:
                            rm_gaps &= set(itertools.chain.from_iterable([list(range(k, entry.gaps[k])) for k in entry.gaps]))
                        rm_gaps = sorted(list(rm_gaps))

                        # make intervals of consecutive gap positions for faster join()
                        gap_ranges = []
                        for k, g in itertools.groupby(enumerate(rm_gaps), lambda x:x[0]-x[1]):
                            group = list(map(itemgetter(1), g))
                            gap_ranges.append((group[0], group[-1])) # tuples with intervals

                        if len(gap_ranges) > 0:
                            if gap_ranges[0][0] != 0:
                                gap_ranges = [(-1, -1)] + gap_ranges

                            if gap_ranges[-1] != lcb.length:
                                gap_ranges = gap_ranges + [(lcb.length, lcb.length)]

                            for entry in entries:
                                entry.sequence = ''.join(
                                    [entry.sequence[(gap_ranges[i][1] + 1):gap_ranges[i + 1][0]] for i in
                                     range(len(gap_ranges) - 1)])
                                if entry.genome_nr > rm_genome:
                                    entry.genome_nr -= 1
                        # if no gaps found only reduce genome nr (avoid looping through entries twice if gaps present)
                        else:
                            for entry in entries:
                                if entry.genome_nr > rm_genome:
                                    entry.genome_nr -= 1

                    elif len(entries) == 1: # if only one entry left replace all gaps in sequence
                        entries[0].sequence = entries[0].sequence.replace("-", "")
                        if entries[0].genome_nr > rm_genome:
                                entries[0].genome_nr -= 1
                else:
                    for entry in entries:
                        if entry.genome_nr > rm_genome:
                            entry.genome_nr -= 1

                lcb.entries[:] = entries

            alignment.lcbs[:] = [lcb for lcb in alignment.lcbs if len(lcb.entries) > 0]

            max_genome = len(alignment.genomes)
            for nr in range(rm_genome + 1, max_genome + 1):
                alignment.genomes[nr - 1] = alignment.genomes[nr]
            del alignment.genomes[max_genome]

            return alignment

        else:
            raise ParameterError("remove_genome", rm_genome,
                                 "between 0 and " + str(len(alignment.genomes)) + " (number of genomes in XMFA)")

    def merge(self, alignment):
        for lcb in alignment.lcbs:
            lcb.entries = sorted(lcb.entries, key=lambda entry: entry.genome_nr)
        new_alignment = Alignment(alignment.xmfa_file)
        for nr, genome in alignment.genomes.items():
            new_alignment.add_genome(genome, nr)
        for order in alignment.genomes:
            if len(alignment.lcbs) > 1:  # if more than one lcb is unchecked try to merge
                alignment.lcbs = alignment.get_sorted_lcbs(order)
                j = 0
                if alignment.lcbs[0].get_entry(order) is not None:
                    new_alignment.add_lcb(alignment.lcbs[0])
                    for lcb in range(1, len(alignment.lcbs)):
                        j += 1
                        if alignment.lcbs[lcb].get_entry(order) is not None:
                            i = 0
                            if len(alignment.lcbs[lcb].entries) == len(alignment.lcbs[lcb - 1].entries):
                                strand = alignment.lcbs[lcb].entries[0].strand == alignment.lcbs[lcb - 1].entries[0].strand
                                for entry in range(0, len(alignment.lcbs[lcb].entries)):
                                    if alignment.lcbs[lcb].entries[entry].genome_nr != alignment.lcbs[lcb - 1].entries[entry].genome_nr \
                                            or alignment.lcbs[lcb].entries[entry].start - alignment.lcbs[lcb - 1].entries[entry].end != 1 \
                                            or (alignment.lcbs[lcb].entries[entry].strand ==alignment.lcbs[lcb - 1].entries[entry].strand) != strand:
                                        # if an entry does not fulfill all conditions stop and do not merge this lcb
                                        new_alignment.add_lcb(alignment.lcbs[lcb])
                                        break
                                    else:
                                        i += 1
                                if i == len(alignment.lcbs[lcb].entries):
                                    if not strand:  # if all entries have an unequal strand reverse complement this lcb
                                        alignment.lcbs[lcb].reverse_complement_entries()
                                    new_alignment.lcbs[-1].length += alignment.lcbs[lcb].length
                                    for pos in range(0, len(new_alignment.lcbs[-1].entries)):
                                        new_alignment.lcbs[-1].entries[pos].sequence += alignment.lcbs[lcb].entries[pos].sequence
                                        new_alignment.lcbs[-1].entries[pos].end = alignment.lcbs[lcb].entries[pos].end
                            else:
                                new_alignment.add_lcb(alignment.lcbs[lcb])
                        else:
                            break
                alignment.lcbs[:] = alignment.lcbs[j:len(alignment.lcbs)]  # continue with unchecked lcbs
            # if only one lcb left check whether it is already checked or not
            elif len(alignment.lcbs) == 1 and new_alignment.lcbs[-1].entries[0].genome_nr != alignment.lcbs[0].entries[0].genome_nr:
                new_alignment.add_lcb(alignment.lcbs[0])  # if not checked yet add it as last lcb and finish
                break
            else:
                break
        return new_alignment