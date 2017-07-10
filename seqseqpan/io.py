import os
import re
import collections
import subprocess
import pdb

import contextlib
import warnings
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

from seqseqpan.exception import *
from seqseqpan.base import *
from seqseqpan.formatter import Splitter


class Parser:
    def parse_xmfa(self, filename):
        alignment = Alignment(filename)

        with open(filename, "r") as xmfa:
            line = xmfa.readline()
            seq_parts = []
            start = 0
            end = 0
            seq_nr = 0
            ses = []
            while line:
                line = line.rstrip()
                if line.startswith("#"):  # parse comment section
                    # each sequence has an associated file (can be the same for more sequences -> multifasta)
                    m = re.match("#Sequence(\d+)File\s+(.+)", line)
                    if m is not None:
                        number = m.group(1)
                        fn = m.group(2)
                        entry = -1
                        line = xmfa.readline()
                        # with multifasta files sequence - entry numbers are reported in line after filename
                        m = re.match("#Sequence" + number + "Entry\s+(\d+)", line)
                        if m is not None:
                            entry = m.group(1)
                            line = xmfa.readline()

                        m = re.match("#Sequence" + number + "Format\s+(\w+)", line)
                        if m is not None:
                            file_format = m.group(1)
                            line = xmfa.readline()
                        genome = Genome(fn, file_format, entry)
                        alignment.add_genome(genome, number)
                        continue

                elif line.startswith(">"):  # a sequence start was encountered
                    if len(seq_parts) > 0 and int(end) > 0:  # save previous sequence
                        seq = "".join(seq_parts)
                        ses.append(SequenceEntry(seq_nr, start, end, strand, seq))
                        alignment.genomes[int(seq_nr)].length = max(alignment.genomes[int(seq_nr)].length, int(end))

                    seq_parts = []
                    m = re.match(">\s*(\d+):(\d+)-(\d+) ([+-]) ", line)
                    if m is not None:
                        seq_nr = m.group(1)
                        start = m.group(2)
                        end = m.group(3)
                        strand = m.group(4)
                    else:
                        raise XMFAHeaderFormatError(line.strip())
                elif line.startswith("="):
                    if len(seq_parts) > 0 and int(end) > 0:
                        seq = "".join(seq_parts)
                        ses.append(SequenceEntry(seq_nr, start, end, strand, seq))
                        alignment.genomes[int(seq_nr)].length = max(alignment.genomes[int(seq_nr)].length, int(end))

                    alignment.add_lcb_entries(ses)

                    seq_parts = []
                    ses = []
                else:
                    seq_parts.append(line)

                line = xmfa.readline()

        return alignment

    def _parse_consensus(self, filename):
        try:
            record = SeqIO.read(open(filename), "fasta")
        except ValueError:
            raise ConsensusFastaError()

        m = re.match("^[^;]+;(\d+)\|(.*)", record.id)

        if m is not None:
            order = m.group(1)
            xmfa_file = m.group(2)
        else:
            raise ConsensusFastaFormatError()

        try:
            cons = Consensus(str(record.seq), order, xmfa_file, filename)
        except ParameterError:
            raise ConsensusFastaFormatError()

        return cons

    def parse_block_separated_consensus(self, filename):
        consensus = self._parse_consensus(filename + ".blockseparated.fasta")
        consensus.block_start_indices = self._parse_consensus_separator(filename)

        return consensus

    def _parse_consensus_separator(self, filename):
        with open(filename + ".blockseparated.idx", "r") as in_file:
            # skip first two lines
            line = in_file.readline()
            line = in_file.readline()
            line = in_file.readline()

            if line != "":
                    return [int(idx) for idx in line.strip().split(";")]
            else:
                    return []

    def parse_consensus_index(self, filename):
        with open(filename + ".idx", "r") as in_file:
            line = in_file.readline()
            m = re.match("#Fasta\t(.+)", line)
            if m is not None:
                fasta_file = m.group(1)
            else:
                raise ConsensusFastaIdxFormatError("Wrong format of Fasta header line.")
            line = in_file.readline()
            m = re.match("#XMFA\t(.+)", line)
            if m is not None:
                xmfa_file = m.group(1)
            else:
                raise ConsensusFastaIdxFormatError("Wrong format of XMFA header line.")

            alignment = Alignment(xmfa_file)

            lcb = LCB()
            lcb_length = 0
            lcb_ends_list = [0]

            line = in_file.readline()
            while line:
                line = line.strip()

                if line.startswith("#"):
                    # each sequence has an associated file (can be the same for more sequences -> multifasta)
                    m = re.match("#Sequence(\d+)File\s+(.+)", line)
                    if m is not None:
                        number = m.group(1)
                        fn = m.group(2)
                        entry = -1
                        line = in_file.readline()
                        # with multifasta files sequence - entry numbers are reported in line after filename
                        m = re.match("#Sequence" + number + "Entry\s+(\d+)", line)
                        if m is not None:
                            entry = m.group(1)
                            line = in_file.readline()

                        m = re.match("#Sequence" + number + "Format\s+(\w+)", line)
                        if m is not None:
                            file_format = m.group(1)
                            line = in_file.readline()
                        genome = Genome(fn, file_format, entry)
                        alignment.add_genome(genome, number)
                        continue

                elif not line == "":

                    fields = line.split("\t")
                    block_id = fields[0]
                    nr = int(fields[1])
                    start = int(fields[2])
                    end = int(fields[3])
                    strand = fields[4]
                    if block_id == "b":
                        # add previous lcb with all entries
                        if lcb_length > 0:
                            lcb.length = lcb_length
                            alignment.add_lcb(lcb)
                            lcb = LCB()

                        # store end of current lcb
                        lcb_length = (end - start) + 1
                        lcb_ends_list.append(end)
                    elif block_id == "s":
                        # add entry to current lcb
                        e = SequenceEntry(nr, start, end, strand, '')

                        if len(fields) == 6:
                            gaps = fields[5]
                            gap_dict = {}
                            for interval in gaps.split(";"):
                                start, end = interval.split("-")
                                gap_dict[int(start)] = int(end)
                            e.gaps = gap_dict
                        lcb.entries.append(e)
                    else:
                        raise ConsensusFastaIdxFormatError("Lines can only start with 'b' or 's'.")
                line = in_file.readline()

            # add last LCB to alignment
            if lcb_length > 0:
                lcb.length = lcb_length
                alignment.add_lcb(lcb)

            consensus = Consensus(sequence="", order=0, xmfa_file=xmfa_file, fasta_file=fasta_file)
            consensus.block_start_indices = lcb_ends_list

        return alignment, consensus

    def parse_mapping_coordinates(self, coord_f):
        with open(coord_f, "r") as in_file:
            header = in_file.readline().strip()
            source, dest = header.split("\t")
            dests = dest.split(",")
            coords = [int(line.strip()) for line in in_file if line.strip() != '']

            if source == "" or dests == "" or len(dests) == 0 or len(coords) == 0:
                raise CoordinatesInputError()

        return source, dests, coords

    def parse_genome_description(self, genome_desc_f):
        chromosome_desc = collections.defaultdict(list)
        with open(genome_desc_f, "r") as in_file:
            line = in_file.readline().strip()
            while line:
                fields = line.split("\t")
                length = int(fields[2]) if len(fields) > 2 else None
                chromosome_desc[int(fields[0])].append({"desc": fields[1], "length": length})
                line = in_file.readline().strip()

        return chromosome_desc

class Writer:
    _mauve_format_string = "#FormatVersion Mauve1\n"
    _mauve_genome_file = '#Sequence{0}File\t{1}\n'
    _mauve_genome_entry = '#Sequence{0}Entry\t{1}\n'
    _mauve_genome_format = '#Sequence{0}Format\t{1}\n'
    _mauve_block_header = '> {0}:{1}-{2} {3} {4}\n'

    _consensus_index_fasta = '#Fasta\t{0}\n'
    _consensus_index_xmfa = '#XMFA\t{0}\n'
    _consensus_index_block_line = '\nb\t{0}\t{1}\t{2}\t+\n'
    _consensus_index_sequence_line = 's\t{0}\t{1}\t{2}\t{3}\t{4}\n'

    _maf_format_string = "##maf version=1\n"
    _maf_sequence_header = "\na label={0}\n"
    _maf_entry_header = "s\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"

    def write_xmfa(self, alignment, path, name, order=0, check_invalid=True):

        if check_invalid and alignment.is_invalid():
            print("\n!!!!!!!!!\n!!!!!!!\nWARNING!!!!!!: XMFA is invalid!\n!!!!!!!!!\n!!!!!!!\n")

        with open(path + "/" + name + ".xmfa", "w") as output:
            output.write(self._mauve_format_string)

            for nr, genome in sorted(alignment.genomes.items()):
                output.write(self._mauve_genome_file.format(nr, genome.file_path))
                if genome.entry > 0:
                    output.write(self._mauve_genome_entry.format(nr, genome.entry))
                output.write(self._mauve_genome_format.format(nr, genome.format))

            sorted_lcbs = alignment.get_sorted_lcbs(order)
            count = 0
            for lcb in sorted_lcbs:
                count += 1
                for entry in sorted(lcb.entries, key=lambda e: e.genome_nr):
                    output.write(self._mauve_block_header.format(entry.genome_nr,
                                                                 entry.start,
                                                                 entry.end,
                                                                 entry.strand,
                                                                 alignment.genomes[entry.genome_nr].file_path
                                                                 )
                                 )
                    output.write("\n".join(re.findall(".{1,80}", entry.sequence)) + "\n")
                output.write("=\n")

    def write_maf(self, alignment, path, name, chromosome_desc, check_invalid=True):

        if check_invalid and alignment.is_invalid():
            print("\n!!!!!!!!!\n!!!!!!!\nWARNING!!!!!!: MAF is invalid!\n!!!!!!!!!\n!!!!!!!\n")

        with open(path + "/" + name + ".maf", "w") as output:
            output.write(self._maf_format_string)

            splitter = Splitter(alignment, chromosome_desc)
            split = splitter.split_alignment()

            count = 0
            for lcb in split.lcbs:

                count += 1
                output.write(self._maf_sequence_header.format(count))

                for entry in sorted(lcb.entries, key=lambda e: e.genome_nr):
                    genome = alignment.genomes[entry.genome_nr]
                    chrom_starts = splitter.get_chromosomes_for_entry(entry)

                    if len(chrom_starts) != 1:
                        raise Exception("Splitting by chromosomes went wrong.")

                    chrom_start = chrom_starts[0]
                    chrom = genome.chromosomes[chrom_start]

                    start = entry.start - chrom_start

                    output.write(self._maf_entry_header.format(chrom["desc"].replace(" ", "_"), start,
                                                               ((entry.end - entry.start) + 1), entry.strand,
                                                               chrom["length"], entry.sequence))

    def write_mapping_coordinates(self, source, destinations, coords_dict, path, name):
        with open(os.path.abspath(path + "/" + name + ".txt"), "w") as output:
            output.write(''.join([str(source), " (source)", "\t"]))
            output.write('\t'.join(destinations))
            output.write("\n")

            for coord, cur_dict in sorted(coords_dict.items()):
                output.write(str(cur_dict[source]) + "\t")
                new_coords = [str(cur_dict.get(dest, "-")) for dest in destinations]
                output.write('\t'.join(new_coords))
                output.write("\n")

    def write_consensus(self, alignment, path, name, order=0):
        filename = os.path.abspath(path + "/" + name + "_consensus.fasta")

        consensus = Consensus()
        consensus.from_alignment(alignment, filename, order)

        self._write_consensus_index(alignment, filename, order)
        self._write_consensus_separator(consensus, alignment, filename)
        header = consensus.get_fasta_header(name)

        record = SeqRecord(Seq(consensus.sequence), id=header, description='')

        with open(filename + ".blockseparated.fasta", "w") as handle:
            SeqIO.write(record, handle, "fasta")

        record.seq = Seq(consensus.get_undelimited_sequence())

        with open(filename, "w") as handle:
            SeqIO.write(record, handle, "fasta")

    def write_fasta(self, seq_name, sequence, path, name):
        filename = os.path.abspath(path + "/" + name + ".fasta")
        record = SeqRecord(Seq(sequence), id=seq_name, description='')

        with open(filename, "w") as handle:
            SeqIO.write(record, handle, "fasta")

    def _write_consensus_index(self, alignment, fastafile, order=0):
        with open(fastafile + ".idx", "w") as output:
            output.write(self._consensus_index_fasta.format(fastafile))
            output.write(self._consensus_index_xmfa.format(alignment.xmfa_file))

            for nr, genome in sorted(alignment.genomes.items()):
                output.write(self._mauve_genome_file.format(nr, genome.file_path))
                if genome.entry > 0:
                    output.write(self._mauve_genome_entry.format(nr, genome.entry))
                output.write(self._mauve_genome_format.format(nr, genome.format))

            sorted_lcbs = alignment.get_sorted_lcbs(order)

            consensus_end = 0
            counter = 1

            for lcb in sorted_lcbs:
                output.write(self._consensus_index_block_line.format(counter, consensus_end + 1,
                                                                     consensus_end + lcb.length))
                consensus_end += lcb.length

                counter += 1
                for entry in lcb.entries:
                    output.write(self._consensus_index_sequence_line.format
                                 (entry.genome_nr, entry.start, entry.end, entry.strand,
                                  ';'.join(
                                      ['-'.join([str(start), str(end)]) for start, end in sorted(entry.gaps.items())])
                                  )
                                 )

    def _write_consensus_separator(self, consensus, alignment, fasta_file):
        with open(fasta_file + ".blockseparated.idx", "w") as output:
            output.write(self._consensus_index_fasta.format(fasta_file))
            output.write(self._consensus_index_xmfa.format(alignment.xmfa_file))
            output.write(';'.join([str(idx) for idx in consensus.block_start_indices]))


class Processor:
    def __init__(self, path, blat="blat"):
        self.blat = blat
        self.path = path

    def external_blat(self, seq_one, seq_two):
        filename_one = os.path.abspath(self.path + "/" + "realigner_realign_seq1.fasta")
        record_one = SeqRecord(Seq(seq_one), id="seq1", description='')

        with open(filename_one, "w") as handle:
            SeqIO.write(record_one, handle, "fasta")

        filename_two = os.path.abspath(self.path + "/" + "realigner_realign_seq2.fasta")
        record_two = SeqRecord(Seq(seq_two), id="seq2", description='')

        with open(filename_two, "w") as handle:
            SeqIO.write(record_two, handle, "fasta")

        output_file = os.path.abspath(self.path + "/" + "realigner_blat_realign.pslx")

        returncode = subprocess.call(args=[self.blat, filename_one, filename_two, output_file, "-out=pslx"])

        try:
            if returncode == 0:
                with open(output_file, "r") as res:
                    line = res.readline()
                    while not line.startswith("-------"):
                        line = res.readline()
                    line = res.readline()
                    if line == '' or line.split(sep='\t')[8] == "-": # query aligned on "-" strand -> cannot be used as realignment!
                        raise ValueError()
                return SearchIO.read(output_file, "blat-psl", pslx=True)
            else:
                raise ValueError()
        except ValueError as e:
            # no alignment --> return None, so full sequences will be used
            return None
        finally:
            with contextlib.suppress(FileNotFoundError):  # in case files were never created...
                pass
                #os.remove(filename_one)
                #os.remove(filename_two)
                #os.remove(output_file)
