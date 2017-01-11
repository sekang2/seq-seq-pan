import sys
import os
import re
import random

from Bio import SeqIO
from Bio.Seq import Seq

from seqseqpan.exception import *
from seqseqpan.constants import BLOCK_DELIMITER, RANDOM_SEED


class Genome:
    def __init__(self, file_path, fileformat, entry=-1):
        self.file_path = os.path.abspath(file_path)
        self.format = fileformat
        self.entry = int(entry)
        self.chromosomes = None
        self.length = 0

    def add_chromosomes(self, chromosomes):
        self.chromosomes = {}
        start = 1
        for chromosome in chromosomes:
            if chromosome["length"] is None:
                chromosome["length"] = self.length
            self.chromosomes[int(start)] = chromosome
            start += chromosome["length"]

class Alignment:
    def __init__(self, xmfa_file):
        self.xmfa_file = os.path.abspath(xmfa_file)
        self.lcbs = []
        self.genomes = {}

    def add_genome(self, genome, number):
        self.genomes[int(number)] = genome

    def add_lcb_entries(self, entries):
        number = len(self.lcbs) + 1
        lcb = LCB(number)
        lcb.add_entries(entries)
        self.lcbs.append(lcb)

    def add_lcb(self, lcb):
        number = len(self.lcbs) + 1
        lcb.number = number
        self.lcbs.append(lcb)

    def get_consensus(self, order=0):
        sorted_lcbs = self.get_sorted_lcbs(order)
        delim = BLOCK_DELIMITER

        consensus_seqs = [lcb.consensus_sequence() for lcb in sorted_lcbs]

        # store beginning of each block calculated using lengths of joined lcbs
        lcb_lengths = [lcb.length for lcb in sorted_lcbs]
        lcb_lengths = [0] + lcb_lengths
        for i in range(1, len(lcb_lengths)):
            lcb_lengths[i] = lcb_lengths[i - 1] + lcb_lengths[i] + len(delim)

        consseq = delim.join(consensus_seqs)

        return consseq, lcb_lengths[:-1]

    def get_sorted_lcbs(self, order):
        # if order == 0 sort by blocknr
        if len(self.genomes) >= order > -1:
            sorted_lcbs = sorted(self.lcbs, key=lambda lcb: lcb.number)

            if order > 0:
                sorted_lcbs = sorted(sorted_lcbs, key=lambda lcb, order=order:
                                     lcb.get_entry(order).start if lcb.get_entry(order) is not None else sys.maxsize)
            return sorted_lcbs
        else:
            raise ParameterError("order", order,
                                 "between 0 and " + str(len(self.genomes)) + " (number of genomes in XMFA)")

    def is_invalid(self):
        invalid = True
        for genome_nr in self.genomes.keys():
            ordered = self.get_sorted_lcbs(genome_nr)
            last_end = 0
            invalid = False
            for lcb in ordered:
                entry = lcb.get_entry(genome_nr)
                if entry is None:
                    break

                invalid = (entry.start - last_end != 1)
                if invalid:
                    break

                last_end = entry.end

            if invalid:
                break
        return invalid


class LCB:
    _a = 1
    _c = 2
    _g = 3
    _t = 4
    _gap = 0
    _rev_table = ("-", "A", "C", "G", "T")
    _iupac_dict = {
        "A": {"value": 4, "pos": [_a]},
        "C": {"value": 4, "pos": [_c]},
        "G": {"value": 4, "pos": [_g]},
        "T": {"value": 4, "pos": [_t]},
        "B": {"value": 2, "pos": [_c, _g, _t]},
        "D": {"value": 2, "pos": [_a, _g, _t]},
        "H": {"value": 2, "pos": [_a, _c, _t]},
        "K": {"value": 3, "pos": [_g, _t]},
        "M": {"value": 3, "pos": [_a, _c]},
        "R": {"value": 3, "pos": [_a, _g]},
        "S": {"value": 3, "pos": [_c, _g]},
        "V": {"value": 2, "pos": [_a, _c, _g]},
        "W": {"value": 3, "pos": [_a, _t]},
        "Y": {"value": 3, "pos": [_c, _t]},
        "-": {"value": 0, "pos": [_gap]},
        "N": {"value": 1, "pos": [_a, _c, _g, _t]}
    }

    def __init__(self, number=0):
        self.number = number
        self.entries = []
        self.length = 0
        random.seed(a=RANDOM_SEED)

    def __eq__(self, other):
        return sorted(self.entries) == sorted(other.entries)

    def add_entries(self, entries):
        if type(entries) is not list:
            entries = [entries]

        length = len(entries[0].sequence)

        # check if new entries have same length as already added ones 
        # and whether all new entries have the same length
        if (self.length > 0 and length != self.length) or (not all(len(e.sequence) == length for e in entries)):
            raise LcbInputError(self.number)

        self.length = length
        self.entries.extend(entries)

    def get_entry(self, genome_nr):
        entry = None
        for e in self.entries:
            if e.genome_nr == genome_nr:
                entry = e
                break
        return entry

    def reverse_complement_entries(self):
        for e in self.entries:
            e.reverse_complement()

    def consensus_sequence(self):
        if len(self.entries) > 0:

            consseq = ''.join([self._majority((lambda i=i: [e.sequence[i] for e in self.entries])())
                               for i in range(len(self.entries[0].sequence))
                               ])  # get i'th letter in sequence of each entry of LCB
            return consseq
        else:
            return ""

    def _majority(self, bases):
        # steps:
        # convert characters into unambiguous code (A,C,G,T)
        # add factor of current character to all positions of unambigious bases
        # 
        # get max position
        # if max is gap and other value greater than zero take that base
        # if value is same as current max change base based on random desicion
        cur_arr = [0, 0, 0, 0, 0]
        try:
            for base in bases:
                base_info = self._iupac_dict[base.upper()]
                for pos in base_info["pos"]:
                    cur_arr[pos] += base_info["value"]
        except KeyError as e:
            raise ConsensusFastaInputError(e.args[0])
        else:
            max_val = 0
            max_idx = 0
            equal = 0
            for i in range(len(cur_arr)):
                cur_val = cur_arr[i]

                if cur_val > max_val or (cur_val > 0 and max_idx == 0):
                    max_val = cur_val
                    max_idx = i
                    equal = 0
                elif cur_val == max_val:
                    equal += 1
                    if random.choice([True, False]):
                        max_idx = i

            if max_val == 0 or equal == (len(cur_arr) - 2):
                return "N"
            else:
                return self._rev_table[max_idx]


class SequenceEntry:
    def __init__(self, genome_nr, start, end, strand, sequence):
        self.genome_nr = int(genome_nr)
        self.sequence = sequence
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __setattr__(self, name, value):
        self.__dict__[name] = value
        if name == "sequence":
            self._gap_dict()

    def __eq__(self, other):
        return (self.sequence == other.sequence and
                self.start == other.start and
                self.strand == other.strand)

    def _gap_dict(self):
        self.gaps = {i.start(): i.end() for i in re.finditer('-+', self.sequence)}

    def get_gap_sublist(self, start=0, end=None):

        if end is None:
            end = len(self.sequence)

        sub_gaps = {gstart: gend for gstart, gend in self.gaps.items() if gend > start and gstart <= end}
        if len(sub_gaps) > 0:
            sorted_keys = sorted(sub_gaps.keys())
            smallest_start = sorted_keys[0]
            highest_start = sorted_keys[-1]

            # if a gap streches over borders of region, make region borders new start or end of gap
            if sub_gaps[highest_start] > end:
                sub_gaps[highest_start] = end + 1
            if smallest_start < start:
                sub_gaps[start] = sub_gaps[smallest_start]
                del sub_gaps[smallest_start]

        return sub_gaps

    def reverse_complement(self):
        self.strand = ("+" if self.strand == "-" else "-")
        seq = Seq(self.sequence)
        self.sequence = str(seq.reverse_complement())

    def get_position_within_entry_with_gaps(self, pos_within_block_without_gaps):
        if pos_within_block_without_gaps == 1:
            return pos_within_block_without_gaps

        cur_nr_nongaps = 0
        pos_within_block = 0

        # go through gaps and add up non-gap characters 
        # stop if wanted position is located before next gap
        # calculate final position

        for start, end in sorted(self.gaps.items()):
            new_nr_nongaps = cur_nr_nongaps + (start - pos_within_block)
            if new_nr_nongaps >= pos_within_block_without_gaps:
                break
            else:
                pos_within_block = end
                cur_nr_nongaps = new_nr_nongaps

        pos_within_block += (pos_within_block_without_gaps - cur_nr_nongaps)

        return pos_within_block


class Consensus:
    def __init__(self, sequence="", order=0, xmfa_file="", fasta_file=""):
        self.sequence = sequence
        self.order = int(order)
        self.xmfa_file = os.path.abspath(xmfa_file)
        self.fasta_file = os.path.abspath(fasta_file)
        self.block_start_indices = []

    def from_alignment(self, alignment, order, fasta_file):
        self.order = int(order)
        self.xmfa_file = alignment.xmfa_file
        self.fasta_file = os.path.abspath(fasta_file)
        self.sequence, self.block_start_indices = alignment.get_consensus(order)

    def get_fasta_header(self, name):
        header = name + ";" + str(self.order) + "|" + self.xmfa_file
        return header

    def get_undelimited_sequence(self):
        seq = self.sequence
        parts = [seq[i:j] for i, j in zip(self.block_start_indices, self.block_start_indices[1:] + [None])]
        delim_len = len(BLOCK_DELIMITER)

        seq = ''.join([part[:-delim_len] for part in parts[:-1]] + [parts[-1]])

        return seq
