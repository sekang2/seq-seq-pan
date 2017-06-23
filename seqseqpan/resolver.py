import bisect
import re

from seqseqpan.exception import ConsensusXMFAInputError, ConsensusCorruptError
from seqseqpan.base import *


class Resolver:
    def resolve_multialignment(self, alignment, consensus, consensus_genome_nr, new_genome_nr):

        if len(alignment.genomes) > 2:
            raise ConsensusXMFAInputError()

        resolved = Alignment(alignment.xmfa_file)
        for nr, genome in alignment.genomes.items():
            resolved.add_genome(genome, nr)

        for lcb in alignment.lcbs:

            consensus_entry = lcb.get_entry(consensus_genome_nr)
            if consensus_entry is not None:
                # consensus sequence entry of LCB should be on forward strand for easier calculation of coordinates
                # also: delimiter positions are stored in forward direction
                # if not: reverse complement all entries in LCB

                if consensus_entry.strand == "-":
                    lcb.reverse_complement_entries()

                # check if there are delimiters in the sequence 
                split_lcbs = self._split_lcb(lcb, consensus_genome_nr, new_genome_nr, consensus)  # split LCB
                for split_lcb in split_lcbs:
                    resolved.add_lcb(split_lcb)
            else:
                resolved.add_lcb(lcb)

        return resolved

    def reconstruct_alignment(self, resolved_align, consensus, orig_align, consensus_genome_nr, new_genome_nr):
        if len(resolved_align.genomes) > 2:
            raise ConsensusXMFAInputError()

        # add genomes from consensus and new genome
        recalculated = Alignment(orig_align.xmfa_file)
        for nr, genome in orig_align.genomes.items():
            recalculated.add_genome(genome, nr)

        added_genome_nr = len(orig_align.genomes) + 1

        # if first genome missing in alignment (singleton alignment!) add new genome to beginning of genome set
        if not (1 in orig_align.genomes):
            added_genome_nr = 1

        recalculated.add_genome(resolved_align.genomes[new_genome_nr], added_genome_nr)
        try:
            sorted_orig_lcbs = orig_align.get_sorted_lcbs(consensus.order)
        except ParameterError:
            raise ConsensusFastaFormatError()
        else:

            # for each block reconstruct original sequence from consensus entry and/or add entry from new genome
            for lcb in resolved_align.lcbs:

                recalculated_lcb = LCB(lcb.number)
                new_entry = lcb.get_entry(new_genome_nr)

                if new_entry is not None:
                    new_entry.genome_nr = added_genome_nr
                    recalculated_lcb.add_entries(new_entry)

                consensus_entry = lcb.get_entry(consensus_genome_nr)
                if consensus_entry is not None:
                    org_entries = self._calculate_coordinates(consensus_entry, consensus, sorted_orig_lcbs)
                    try:
                        recalculated_lcb.add_entries(org_entries)
                    except LcbInputError as e:
                        e.message += " (Error occurred in recalculating coordinates step.)"
                        raise e
                    except IndexError as e:
                        print(lcb.number)
                        raise e

                recalculated.add_lcb(recalculated_lcb)

            return recalculated

    def _split_lcb(self, lcb, consensus_genome_nr, new_genome_nr, consensus):
        split_lcbs = []
        delim_positions = []

        consensus_entry = lcb.get_entry(consensus_genome_nr)

        # check whether consensus_entry overlaps delimiter position

        # check if entry contains delimiter sequence less than normal length at beginning
        intervals = [idx - 1000 for idx in sorted(consensus.block_start_indices)[1:] if
                     0 <= (idx - consensus_entry.start) <= len(BLOCK_DELIMITER)]

        # check for all overlaps
        intervals.extend([idx - 1000 for idx in sorted(consensus.block_start_indices)[1:] if
                          consensus_entry.start < (idx - 1000) < consensus_entry.end])

        for startInterval in intervals:
            start_within_block = startInterval - consensus_entry.start + 1
            offset = (start_within_block * -1 if start_within_block < 0 else 0)
            # in case of delimiter sequencing at beginning of entry, index will be negative but should be 0
            start_within_block = (consensus_entry.get_position_within_entry_with_gaps(start_within_block)
                                  if start_within_block > 0 else 0)

            # search for delimiter sequence starting at interval position:
            # N with gaps inbetween with maximum length of 1000
            m = re.search('(N-*){0,' + str(len(BLOCK_DELIMITER) - offset - 1) + '}N',
                          consensus_entry.sequence[start_within_block:])
            if m is None:
                raise ConsensusCorruptError(start_within_block + consensus_entry.start)

            # store positions to split
            delim_positions.append(start_within_block + m.start(0))
            delim_positions.append(start_within_block + m.end(0))

        if len(delim_positions) == 0:  # no entries with delimiter sequence in LCB
            return [lcb]

        indices = delim_positions

        # store whether even or uneven elements are delimiters
        even_is_delimiter = True
        if delim_positions[0] > 0:
            indices = [0] + indices
            even_is_delimiter = False

        if delim_positions[-1] < lcb.length:
            indices = indices + [lcb.length]

        for i in range(len(indices) - 1):

            new_entry = lcb.get_entry(new_genome_nr)

            if new_entry is None:
                entries = []
            else:
                entries = [self._get_split_entry(new_entry, indices[i], indices[i + 1])]

            # check if element is non-delimiter element
            if (i % 2 == 0 and not even_is_delimiter) or (i % 2 != 0 and even_is_delimiter):
                entries.append(self._get_split_entry(consensus_entry, indices[i], indices[i + 1]))

            entries = list(filter(None, entries))  # add all entries that are not None)

            if len(entries) > 0:
                # only one entry - sequence of second genome not present in current LCB --> remove all gaps
                if len(entries) == 1:
                    entries[0].sequence = entries[0].sequence.replace("-", "")
                split_lcb = LCB()
                try:
                    split_lcb.add_entries(entries)
                except LcbInputError as e:
                    raise LcbInputError(e.message + " (Error occured in splitting step.)")
                else:
                    split_lcbs.append(split_lcb)

        return split_lcbs

    def _get_split_entry(self, entry, split_start, split_end):
        seq = entry.sequence[split_start:split_end]
        if seq == "" or re.search("[^-]", seq) is None:
            return None
        split_len = split_end - split_start
        return self._get_new_entry(seq, split_start, split_len, entry)


    def _calculate_coordinates(self, consensus_entry, consensus, orig_lcb_list):

        # calculate in which consensus block this entry is located
        idx = bisect.bisect_left(consensus.block_start_indices, consensus_entry.start)

        idx -= 1
        orig_lcb = orig_lcb_list[idx]

        # calculate start and end of sequence within current consensus block
        start_within_block = consensus_entry.start - consensus.block_start_indices[idx] - 1

        end_within_block = consensus_entry.end - consensus.block_start_indices[idx] - 1
        len_within_block = end_within_block - start_within_block + 1

        # count gaps in consensus entry
        sum_gaps = 0
        if len(consensus_entry.gaps) > 0:
            sum_gaps = sum(end - start for start, end in consensus_entry.gaps.items())

        # get sequence of original blocks the same length as the consensus entry sequence minus number of gaps
        end_sub_sequence = start_within_block + len(consensus_entry.sequence) - sum_gaps

        org_entries = []
        for e in orig_lcb.entries:

            sequence = e.sequence[start_within_block:end_sub_sequence]

            new_entry = self._get_new_entry(sequence,start_within_block,len_within_block,e)

            # add all gaps in consensus sequence to original sequence for correct alignment
            for gap_start, gap_end in sorted(consensus_entry.gaps.items()):
                sequence = self._insert_gap(sequence, gap_start, gap_end - gap_start)

            if sequence != "" and re.search("[^-]", sequence) is not None:
                new_entry.sequence = sequence
                org_entries.append(new_entry)

        return org_entries

    def _insert_gap(self, sequence, position, length):
        seq = ("-" * length).join([sequence[:position], sequence[position:]])
        return seq


    def _get_new_entry(self, split_entry_sequence, split_start, split_len, org_entry):
        entry_len = org_entry.end - org_entry.start + 1

        # count gaps in original sequence before start of sequence to get correct start and end for entry
        sum_gaps_before_split_entry = 0

        subgaps = org_entry.get_gap_sublist(0, split_start - 1)
        if len(subgaps) > 0:
            sum_gaps_before_split_entry = sum(end - start for start, end in subgaps.items())

        if org_entry.strand == "+":
            start = org_entry.start + split_start - sum_gaps_before_split_entry
            end = start + split_len - 1
        else:
            end = org_entry.start + entry_len - split_start + sum_gaps_before_split_entry - 1
            start = end - split_len + 1

        new_entry = SequenceEntry(org_entry.genome_nr, start, end, org_entry.strand, split_entry_sequence)
        sum_gaps_in_split_entry = sum(end - start for start, end in new_entry.gaps.items())
        if org_entry.strand == "+":
            new_entry.end = new_entry.end - sum_gaps_in_split_entry
        else:
            new_entry.start = new_entry.start + sum_gaps_in_split_entry

        return new_entry