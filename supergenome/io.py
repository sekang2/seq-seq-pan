import os
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from supergenome.exception import *
from supergenome.base import *
from supergenome.formatter import Splitter

class Parser:
    def parseXMFA(self, filename):
        alignment = Alignment(filename)
        
        with open(filename, "r") as xmfa:
            line = xmfa.readline()
            seq = ""
            seqparts = []
            start = 0
            end = 0
            seqNr = 0
            ses = []
            while line:
                line = line.rstrip()
                if line.startswith("#"): #parse comment section
                    m = re.match("#Sequence(\d+)File\s+(.+)", line) # each sequence has an associated file (can be the same for more sequences -> multifasta)
                    if m is not None:
                        number = m.group(1)
                        fn = m.group(2)
                        entry = -1
                        line = xmfa.readline()
                        m = re.match("#Sequence"+number+"Entry\s+(\d+)", line) # with multifasta files sequence - entry numbers are reported in line after filename
                        if m is not None:
                            entry = m.group(1)
                            line = xmfa.readline()
                            
                        m = re.match("#Sequence"+number+"Format\s+(\w+)", line) 
                        if m is not None:
                            format = m.group(1)
                            line = xmfa.readline()
                        genome = Genome(fn, format, entry)
                        alignment.addGenome(genome, number)
                        continue
                
                elif line.startswith(">"): # a sequence start was encountered
                    if len(seqparts) > 0 and int(end) > 0: # save previous sequence
                        seq = "".join(seqparts)
                        ses.append(SequenceEntry(seqNr, start, end, strand, seq))
                        
                    seqparts = []
                    m = re.match(">\s*(\d+):(\d+)-(\d+) ([+-]) ", line)
                    if m is not None:
                        seqNr = m.group(1)
                        start = m.group(2)
                        end = m.group(3)
                        strand = m.group(4)
                    else:
                        raise XMFAHeaderFormatError(line.strip())
                elif line.startswith("="):
                    seq = "".join(seqparts)
                    
                    ses.append(SequenceEntry(seqNr, start, end, strand, seq))
                    
                    alignment.addLCBEntries(ses)
                    
                    seqparts = []
                    ses = []
                else:
                    seqparts.append(line)
                    
                line = xmfa.readline()
                    
        return alignment
    
    
    def _parseConsensus(self, filename):
        try:
            record = SeqIO.read(open(filename), "fasta")
        except ValueError as e:
            raise ConsensusFastaError()
        
        m = re.match("^[^;]+;(\d+)\|(.*)", record.id)
        
        if m is not None:
            order = m.group(1)
            xmfaFile = m.group(2)
        else:
            raise ConsensusFastaFormatError()

        try:
            cons = Consensus(str(record.seq), order, xmfaFile, filename)
        except ParameterError:
            raise ConsensusFastaFormatError()
            
        return cons
    
    
    def parseBlockSeparatedConsensus(self, filename):
        consensus = self._parseConsensus(filename+".blockseparated.fasta")
        consensus.blockStartIndices = self._parseConsensusSeparator(filename)
        
        return consensus
    
    
    def _parseConsensusSeparator(self, filename):
        with open(filename+".blockseparated.idx", "r") as input:
            line = input.readline()
            line = input.readline()
            line = input.readline()

            blockStartIndices = [int(idx) for idx in line.strip().split(";")]
            
        return blockStartIndices
    
    
    def parseConsensusIndex(self, filename):
        with open(filename+".idx", "r") as input:
            line = input.readline()
            m = re.match("#Fasta\t(.+)", line)
            if m is not None:
                fastaFile = m.group(1)
            else:
                raise ConsensusFastaIdxFormatError("Wrong format of Fasta header line.")
            line = input.readline()
            m = re.match("#XMFA\t(.+)", line)
            if m is not None:
                xmfaFile = m.group(1)
            else:
                raise ConsensusFastaIdxFormatError("Wrong format of XMFA header line.")
            
            alignment = Alignment(xmfaFile)
            
            lcb = LCB()
            lcbLength = 0
            lcbEndsList = [0]
            
            line = input.readline()
            while line:
                line = line.strip()
                
                if line.startswith("#"):
                    m = re.match("#Sequence(\d+)File\s+(.+)", line) # each sequence has an associated file (can be the same for more sequences -> multifasta)
                    if m is not None:
                        number = m.group(1)
                        fn = m.group(2)
                        entry = -1
                        line = input.readline()
                        m = re.match("#Sequence"+number+"Entry\s+(\d+)", line) # with multifasta files sequence - entry numbers are reported in line after filename
                        if m is not None:
                            entry = m.group(1)
                            line = input.readline()
                            
                        m = re.match("#Sequence"+number+"Format\s+(\w+)", line) 
                        if m is not None:
                            format = m.group(1)
                            line = input.readline()
                        genome = Genome(fn, format, entry)
                        alignment.addGenome(genome, number)
                        continue
                
                elif not line == "": 
                    
                    fields = line.split("\t")
                    id = fields[0]
                    nr = int(fields[1])
                    start = int(fields[2])
                    end = int(fields[3])
                    strand = fields[4]
                    if id == "b":
                        # add previous lcb with all entries
                        if lcbLength > 0:
                            lcb.length = lcbLength
                            alignment.addLCB(lcb)
                            lcb = LCB()
                            
                        # store end of current lcb
                        lcbLength = (end - start) + 1
                        lcbEndsList.append(end)
                    elif id == "s":
                        # add entry to current lcb
                        e = SequenceEntry(nr, start, end, strand, '')
                        
                        if len(fields) == 6:
                            gaps = fields[5]
                            gapDict = {}
                            for interval in gaps.split(";"):
                                start, end = interval.split("-")
                                gapDict[int(start)] =  int(end)
                            e.gaps = gapDict
                        lcb.entries.append(e)
                    else:
                        raise ConsensusFastaIdxFormatError("Lines can only start with 'b' or 's'.")
                line = input.readline()
                
            # add last LCB to alignment
            if lcbLength > 0:
                lcb.length = lcbLength
                alignment.addLCB(lcb)
            
            consensus = Consensus( sequence="", order=0, xmfaFile=xmfaFile, fastaFile=fastaFile)
            consensus.blockStartIndices = lcbEndsList
                
        return (alignment, consensus)
        
        
    def parseMappingCoordinates(self, coord_f):
        with open(coord_f) as input:
            
            header = input.readline().strip()
            source, dest = header.split("\t")
            dests = dest.split(",")
            coords = [int(line.strip()) for line in input]
            
            if source == "" or dests == "" or len(dests) == 0 or len(coords) == 0:
                raise CoordinatesInputError()
        
        return source, dests, coords


        
        
class Writer:
        _mauveFormatString = "#FormatVersion Mauve1\n"
        _mauveGenomeFile = '#Sequence{0}File\t{1}\n'
        _mauveGenomeEntry = '#Sequence{0}Entry\t{1}\n'
        _mauveGenomeFormat = '#Sequence{0}Format\t{1}\n'
        _mauveBlockHeader = '> {0}:{1}-{2} {3} {4}\n'
        
        _consensusIndexFasta = '#Fasta\t{0}\n'
        _consensusIndexXmfa = '#XMFA\t{0}\n'
        _consensusIndexBlockLine = '\nb\t{0}\t{1}\t{2}\t+\n'
        _consensusIndexSequenceLine = 's\t{0}\t{1}\t{2}\t{3}\t{4}\n'
        
        _mafFormatString = "##maf version=1\n"
        _mafSequenceHeader = "\na label={0}\n"
        _mafEntryHeader = "s\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
        
        
        def writeXMFA(self, alignment, path, name, order=0):
            
            if alignment.isInvalid():
                print("\n!!!!!!!!!\n!!!!!!!\nWARNING!!!!!!: XMFA is invalid!\n!!!!!!!!!\n!!!!!!!\n")
                    
        
            with open(path+"/"+name+".xmfa", "w") as output:
                output.write(self._mauveFormatString)
                               
                for nr, genome in sorted(alignment.genomes.items()):
                    output.write(self._mauveGenomeFile.format(nr, genome.filepath))
                    if genome.entry > 0 :
                        output.write(self._mauveGenomeEntry.format(nr, genome.entry))
                    output.write(self._mauveGenomeFormat.format(nr, genome.format))
                
                sortedLCBs = alignment.getSortedLCBs(order)
                count = 0
                for lcb in sortedLCBs:
                    count += 1
                    for entry in sorted(lcb.entries, key=lambda e: e.genomeNr):
                        output.write(self._mauveBlockHeader.format( entry.genomeNr, 
                                                                   entry.start, 
                                                                   entry.end, 
                                                                   entry.strand, 
                                                                   alignment.genomes[entry.genomeNr].filepath
                                                                  )
                                    )
                        output.write("\n".join(re.findall(".{1,80}", entry.sequence))+"\n")
                    output.write("=\n")
        
        
        def writeMAF(self, alignment, path, name, order=0):
            
            if alignment.isInvalid():
                print("\n!!!!!!!!!\n!!!!!!!\nWARNING!!!!!!: MAF is invalid!\n!!!!!!!!!\n!!!!!!!\n")
            
            with open(path+"/"+name+".maf", "w") as output:
                output.write(self._mafFormatString)
                
                sortedLCBs = alignment.getSortedLCBs(order)
                
                splitter = Splitter(alignment)
                splittedLCBs = []
                for lcb in sortedLCBs:
                    splittedLCBs.extend(splitter.splitByChromosomes(lcb))
                
                count = 0
                for lcb in splittedLCBs:
                    
                    count += 1
                    output.write(self._mafSequenceHeader.format(count))
                    
                    for entry in sorted(lcb.entries, key=lambda e: e.genomeNr):
                        genome = alignment.genomes[entry.genomeNr]
                        chrstarts = splitter.getChromosomesForEntry(entry)
                        
                        if len(chrstarts) != 1:
                            raise Exception("Splitting by chromosomes went wrong.")
                        
                        chrstart = chrstarts[0]
                        chr = genome.chromosomes[chrstart]
                        
                        start = entry.start - chrstart
                        
                        output.write(self._mafEntryHeader.format(chr["desc"].replace(" ", "_"), start, ((entry.end - entry.start)+1), entry.strand, chr["length"], entry.sequence))
                    
        
        def writeMappingCoordinates(self, source, dests, coords_dict, path, name):
            with open(os.path.abspath(path + "/"+ name + ".txt"), "w") as output:
                output.write(''.join([str(source)," (source)", "\t"]))
                output.write('\t'.join(dests))
                output.write("\n")
                
                for coord, cur_dict in sorted(coords_dict.items()):
                    output.write(str(coord) + "\t")
                    new_coords = [str(cur_dict.get(dest,"-")) for dest in dests]
                    output.write('\t'.join(new_coords))
                    output.write("\n")

        
        def writeConsensus(self, alignment, unambiguous, path, name, order=0):
            filename = os.path.abspath(path+"/"+name+"_consensus.fasta")
            
            consensus = Consensus()
            consensus.fromAlignment(alignment, order, filename, unambiguous)
            
            self._writeConsensusIndex(alignment, filename, order)
            self._writeConsensusSeparator(consensus, alignment, filename, order)
            header = consensus.getFastaHeader(name)
            
            record = SeqRecord(Seq(consensus.sequence), id=header, description='')
            
            with open(filename + ".blockseparated.fasta", "w") as handle:
                SeqIO.write(record, handle, "fasta")
            
            record.seq = Seq(consensus.getUndelimitedSequence())
            
            with open(filename, "w") as handle:
                SeqIO.write(record, handle, "fasta")
            
        
        def _writeConsensusIndex(self, alignment, fastafile, order=0):
            with open(fastafile+".idx", "w") as output:
                output.write(self._consensusIndexFasta.format(fastafile))
                output.write(self._consensusIndexXmfa.format(alignment.xmfaFile))
                
                for nr, genome in sorted(alignment.genomes.items()):
                    output.write(self._mauveGenomeFile.format(nr, genome.filepath))
                    if genome.entry > 0 :
                        output.write(self._mauveGenomeEntry.format(nr, genome.entry))
                    output.write(self._mauveGenomeFormat.format(nr, genome.format))
                
                sortedLCBs = alignment.getSortedLCBs(order)
            
                consensusEnd = 0
                counter = 1
            
                for lcb in sortedLCBs:
                    output.write(self._consensusIndexBlockLine.format(counter, consensusEnd+1, consensusEnd+lcb.length))
                    consensusEnd += (lcb.length)
                    
                    counter += 1
                    for entry in lcb.entries:
                        output.write(self._consensusIndexSequenceLine.format
                                            ( entry.genomeNr, entry.start, entry.end, entry.strand,
                                              ';'.join(['-'.join( [str(start), str(end)]) for start, end in sorted(entry.gaps.items()) ])
                                            )
                                    )
        
        
        def _writeConsensusSeparator(self, consensus, alignment, fastafile, order):
            with open(fastafile+".blockseparated.idx", "w") as output:
                output.write(self._consensusIndexFasta.format(fastafile))
                output.write(self._consensusIndexXmfa.format(alignment.xmfaFile))
                output.write(';'.join([str(idx) for idx in consensus.blockStartIndices]))
                
