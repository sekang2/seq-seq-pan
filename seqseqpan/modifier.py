import collections
import itertools
from operator import itemgetter
import math

from Bio import pairwise2

from seqseqpan.base import *


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
            use_prev = False
            use_next = False
            prev_gap = False
            next_gap = False
            to_reverse_complement = False

            next_new_entry = None
            prev_new_entry = None

            # check if new entry is small and only created by splitting or aligning of new genome (consensus is None)
            if consensus_entry is None and new_entry is not None and len(new_entry.sequence) <= block_length:

                nr_gaps = len(new_entry.sequence)

                # check if can be appended to previous entry and if that ends with a gap
                if len(merged_split_lcbs) > 0:
                    prev_lcb = merged_split_lcbs[-1]
                    prev_new_entry = prev_lcb.get_entry(new_genome_nr)
                    if prev_new_entry is not None and (new_entry.start - prev_new_entry.end) == 1:
                        use_prev = True
                        if (prev_new_entry.strand == "+" and prev_new_entry.sequence[-1] == "-") \
                                or ( prev_new_entry.strand == "-" and prev_new_entry.sequence[0] == "-"):
                            prev_gap = True


                # check if can be prepended to next entry and if that starts with a gap
                if len(lcbs) > (i+1):
                    next_lcb = lcbs[i + 1]
                    next_new_entry = next_lcb.get_entry(new_genome_nr)
                    if next_new_entry is not None and (next_new_entry.start - new_entry.end) == 1:
                        use_next = True
                        if (next_new_entry.strand == "+" and next_new_entry.sequence[0] == "-")\
                                or (next_new_entry.strand == "-" and next_new_entry.sequence[-1] == "-"):
                            next_gap = True


                neighbour_new_entry = None
                # if both, choose the one with gap at start or end
                if (not next_gap) and use_prev:

                    neighbour_new_entry = prev_new_entry
                    neighbour_lcb = prev_lcb

                    if neighbour_new_entry.strand == "+":
                        append = True
                    else:
                        prepend = True

                elif use_next:
                    neighbour_new_entry = next_new_entry
                    neighbour_lcb = next_lcb
                    if neighbour_new_entry.strand == "+":
                        prepend = True
                    else:
                        append = True

                if neighbour_new_entry is not None:
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
    # local realignment around overlapping or consecutive gaps in two sequences
    # if border_aln_length is not None align only sequences at block border up to given length

    def realign(self, alignment, processor, border_aln_length=0):
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
                    seq_one, seq_two = self._realign(entry_one.sequence, entry_two.sequence, one_first_two_second,
                                                     processor, border_aln_length)
                    entry_one.sequence = seq_one
                    entry_two.sequence = seq_two

                # get regions to realign for updated entries with second entry first
                two_first_one_second = self._get_realign_regions(entry_two.gaps, entry_one.gaps)

                if len(two_first_one_second) > 0:
                    seq_two, seq_one = self._realign(entry_two.sequence, entry_one.sequence, two_first_one_second,
                                                     processor, border_aln_length)
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

    def _realign(self, seq_one, seq_two, realign_regions, processor, border_aln_length):
        # if border_aln_length is not None align only sequences at block border up to given length
        realign_regions = sorted(realign_regions)

        index_offset = 0

        for interval in realign_regions:
            max_index, max_seq_length = max(
                enumerate([interval[0][1] - interval[0][0], interval[1][1] - interval[1][0]]), key=lambda p: p[1])

            # get length of gaps at start or end of block
            short_border_intervals = [(i[1] - i[0]) <= border_aln_length for i in interval  # check border length
                                      if i[0] == 0                                      # start of block
                                      or ((i[1] - index_offset) == len(seq_one))        # end of block
                                      or ((i[1] - index_offset) == len(seq_two))]       # end of block

            interval_start = interval[max_index][0] - index_offset
            interval_end = interval[max_index][1] - index_offset

            # border_aln_length not set OR small sequence at start or end of block
            if border_aln_length == 0 or any(short_border_intervals):

                # check if interval only 'N' - if yes: do not realign
                n_stretch = 'N' * max_seq_length
                if not (seq_one[interval_start:interval_end] == n_stretch or
                                seq_two[interval_start:interval_end] == n_stretch):

                    max_seq_length = math.ceil(max_seq_length * 1.5)

                    # get surrounding sequences
                    seq_start = interval_start - max_seq_length
                    seq_end = interval_end + max_seq_length

                    # do not go over boundaries of sequences!
                    seq_start = max(seq_start, 0)
                    min_orig_seq_length = min(len(seq_one), len(seq_two))
                    seq_end = min(seq_end, min_orig_seq_length)

                    # N-stretches in sequences
                    #  find N-stretch between start and interval and start sub-sequence after nearest stretch
                    n_stretch_length = 10
                    n_stretch = 'N' * n_stretch_length
                    n_stretch_idx = seq_one.rfind(n_stretch, seq_start, interval_start)

                    if n_stretch_idx > -1:
                        seq_start = max(seq_start, (n_stretch_idx + n_stretch_length))
                    n_stretch_idx = seq_two.rfind(n_stretch, seq_start, interval_start)
                    if n_stretch_idx > -1:
                        seq_start = max(seq_start, (n_stretch_idx + n_stretch_length))

                    # find N-stretch between interval and end  and end sub-sequence before nearest stretch
                    n_stretch_idx_one = seq_one.find(n_stretch, interval_end, seq_end)
                    n_stretch_idx_two = seq_two.find(n_stretch, interval_end, seq_end)

                    seq_end = min(seq_end,
                                  (n_stretch_idx_one if n_stretch_idx_one > -1 else seq_end),
                                  (n_stretch_idx_two if n_stretch_idx_two > -1 else seq_end)
                                  )

                    #if border_aln_length > 0:
                    #    print(seq_one[seq_start:seq_end])
                    #    print(seq_two[seq_start:seq_end])

                    seq_one_nogap = seq_one[seq_start:seq_end].replace("-", "")
                    seq_two_nogap = seq_two[seq_start:seq_end].replace("-", "")

                    if not (seq_one_nogap == '' or seq_two_nogap == ''):  # else: do nothing for current interval
                        if (seq_end - seq_start) < 1000: #https://github.com/biopython/biopython/pull/782
                            alignments = pairwise2.align.globalxs(seq_one_nogap.upper(),
                                                                  seq_two_nogap.upper(),
                                                                  -0.5, -0.1, #default
                                                                  one_alignment_only=True)
                            if len(alignments) > 0:
                                max_score = max([x[2] for x in alignments])
                                alignments = (lambda max_score=max_score: [item for item in alignments if item[2] == max_score])()

                                max_length = min([x[4] for x in alignments])
                                alignments = (lambda max_length=max_length: [item for item in alignments if item[4] == max_length])()

                                aln_seq_one = alignments[0][0]
                                aln_seq_two = alignments[0][1]


                            else:
                                # no alignment, do nothing for current interval
                                break

                        else:
                            # do external blat alignment
                            aln_seq_one, aln_seq_two = self._external_blat(seq_one_nogap, seq_two_nogap, processor)


                            if aln_seq_one is not None:
                                max_length = len(aln_seq_one)
                            else:
                                # no alignment, do nothing for current interval
                                break

                        if max_length > 0 and max_length < ( seq_end - seq_start ):  # only use new alignment if better than old one
                            #if border_aln_length > 0:
                            #    print(aln_seq_one)
                            #    print(aln_seq_two)

                            seq_one = aln_seq_one.join([seq_one[:seq_start], seq_one[seq_end:]])
                            seq_two = aln_seq_two.join([seq_two[:seq_start], seq_two[seq_end:]])

                            index_offset += ((seq_end - seq_start) - max_length)

        return seq_one, seq_two

    def _external_blat(self, seq_one, seq_two, processor):

        def _check_tuple_overlap(tuple_list):
            num_tuple = len(tuple_list)

            last_ok_end = -1

            exclude_idx = []

            #overlap_detected = False
            for idx in range(num_tuple):
                if tuple_list[idx][0] < last_ok_end:  #current start smaller than last end -> exclude
                    exclude_idx.append(idx)
                else:
                    last_ok_end = tuple_list[idx][1]  # current tuple is ok -> store end

            return exclude_idx


        align_result = processor.external_blat(seq_one, seq_two)
        if align_result is not None:
            # seq_one equals hit   seq_two equals query

            # only one query and first hit has highest score per definition of blat
            match = align_result[0][0]

            # add starting sequences, in case query or hit do not start at "0"
            seq_one_list = []
            seq_two_list = []

            if match.hit_start > 0:
                seq_one_list.append(seq_one[0:match.hit_start])
                seq_two_list.append("-" * match.hit_start)

            if match.query_start > 0:
                seq_one_list.append("-" * match.query_start)
                seq_two_list.append(seq_two[0:match.query_start])

            if match.is_fragmented:

                # in really rare cases fragments are overlapping!
                exclude_keys = _check_tuple_overlap(match.hit_range_all)
                exclude_keys = exclude_keys + _check_tuple_overlap(match.query_range_all)

                hspfrag_keys = list(set(range(len(match.hit_range_all))) - set(exclude_keys))

                hspfrag_key_num = len(hspfrag_keys)

                for hspfrag_idx in range(hspfrag_key_num):

                    hspfrag_key = hspfrag_keys[hspfrag_idx]
                    hspfrag = match[hspfrag_key]

                    seq_one_list.append(str(hspfrag.hit.seq))
                    seq_two_list.append(str(hspfrag.query.seq))

                    # add sequences between aligned intervals to sequences
                    if hspfrag_idx < (hspfrag_key_num - 1):
                        next_hspfrag = match[hspfrag_keys[hspfrag_idx + 1]]
                        inter_hit_len = next_hspfrag.hit_start - hspfrag.hit_end

                        if inter_hit_len > 0:
                            seq_one_list.append(seq_one[hspfrag.hit_end:next_hspfrag.hit_start])
                            seq_two_list.append("-" * inter_hit_len)

                        inter_query_len = next_hspfrag.query_start - hspfrag.query_end
                        if inter_query_len > 0:
                            seq_one_list.append("-" * inter_query_len)
                            seq_two_list.append(seq_two[hspfrag.query_end:next_hspfrag.query_start])

            else:
                seq_rec_one = match.aln[0]

                if seq_rec_one.id == "seq1":
                    seq_one_list.append(str(match.aln[0].seq))
                    seq_two_list.append(str(match.aln[1].seq))
                else:
                    seq_one_list.append(str(match.aln[1].seq))
                    seq_two_list.append(str(match.aln[0].seq))

            # add last sequence parts if hit or query do not include sequence ends
            if match.hit_end < len(seq_one):
                seq_len = len(seq_one) - match.hit_end
                seq_one_list.append(seq_one[match.hit_end:])
                seq_two_list.append("-" * seq_len)

            if match.query_end < len(seq_two):
                seq_len = len(seq_two) - match.query_end
                seq_one_list.append("-" * seq_len)
                seq_two_list.append(seq_two[match.query_end:])

            return "".join(seq_one_list).upper(), "".join(seq_two_list).upper()
        else:
            return None, None


class SingletonAligner:
    def genome_count_split(self, alignment):
        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()

        single_alignment_1 = Alignment(alignment.xmfa_file)
        single_alignment_2 = Alignment(alignment.xmfa_file)
        single_alignment_1.genomes[1] = alignment.genomes[1]
        single_alignment_2.genomes[2] = alignment.genomes[2]

        pairlcbs = []
        for lcb in alignment.lcbs:
            if len(lcb.entries) == 1:
                if lcb.entries[0].strand == "-":
                    lcb.reverse_complement_entries()

                if lcb.entries[0].genome_nr == 1:
                    single_alignment_1.add_lcb(lcb)
                elif lcb.entries[0].genome_nr == 2:
                    single_alignment_2.add_lcb(lcb)
            else:
                pairlcbs.append(lcb)
        alignment.lcbs = pairlcbs

        return alignment, single_alignment_1, single_alignment_2

    def join(self, alignment, alignment_two):

        if len(alignment.genomes) > 2:
            print("ERROR: I don't want to think about cases with more than 2 genomes now.")

        genome_names_one = [genome.file_path for genome in alignment.genomes.values()]
        genome_names_two = [genome.file_path for genome in alignment_two.genomes.values()]

        if genome_names_one.sort() != genome_names_two.sort():
            print("ERROR: Can only join alignments from same set of genomes.")

        if alignment.genomes[1].file_path != alignment_two.genomes[1].file_path:
            genome_nr_dict = {1: 2, 2: 1}
        else:
            genome_nr_dict = {1: 1, 2: 2}

        for lcb in alignment_two.lcbs:
            for entry in lcb.entries:
                entry.genome_nr = genome_nr_dict[entry.genome_nr]

            alignment.add_lcb(lcb)

        return alignment


class Remover:
    def remove(self, alignment, rm_genome):
        if len(alignment.genomes) >= rm_genome > -1:

            for lcb in alignment.lcbs:

                entries = [entry for entry in lcb.entries if entry.genome_nr != rm_genome]

                # did LCB include entry of genome to remove?
                if len(entries) < len(lcb.entries):
                    # are there any entries left in LCB?
                    if len(entries) > 1:
                        # if more than one entry left search for gaps that are present in all remaining entries

                        rm_gaps = set(itertools.chain.from_iterable(
                            [list(range(k, entries[0].gaps[k])) for k in entries[0].gaps]))

                        for entry in entries[1:]:
                            rm_gaps &= set(
                                itertools.chain.from_iterable([list(range(k, entry.gaps[k])) for k in entry.gaps]))
                        rm_gaps = sorted(list(rm_gaps))

                        # make intervals of consecutive gap positions for faster join()
                        gap_ranges = []
                        for k, g in itertools.groupby(enumerate(rm_gaps), lambda x: x[0] - x[1]):
                            group = list(map(itemgetter(1), g))
                            gap_ranges.append((group[0], group[-1]))  # tuples with intervals

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

                    elif len(entries) == 1:  # if only one entry left replace all gaps in sequence
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
                                            or (alignment.lcbs[lcb].entries[entry].strand == alignment.lcbs[lcb - 1].entries[entry].strand) != strand:

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
