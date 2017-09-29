import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from seqseqpan.io import Parser
from seqseqpan.mapper import Mapper


def main():
    global args
    parser = Parser()
    mapper = Mapper()

# map part

    sparse_align, sparse_consensus = parser.parse_consensus_index(args.consensus_f)
    source, destinations_org, coordinates = parser.parse_mapping_coordinates(args.coord_f)
    destinations = list(set(destinations_org))
    coords_dict = mapper.map_coordinates(sparse_align, sparse_consensus,
                                                         source, destinations_org, coordinates)
# open files part

    print("start open files")
    with open(os.path.abspath(args.output_p + "/" + args.output_name + ".txt"), "w") as output:
        fasta_dict = {}
        if source == 'c':
            fasta_dict[source] = SeqIO.read(open(args.consensus_f), 'fasta')
        else:
            fasta_dict[source] = SeqIO.read(open(sparse_align.genomes[int(source)].file_path), 'fasta')
        for genome in destinations:
            if not genome == 'c':
                fasta_dict[genome] = SeqIO.read(open(sparse_align.genomes[int(genome)].file_path), 'fasta')
            else:
                fasta_dict[genome] = SeqIO.read(open(args.consensus_f), 'fasta')
        print("end open files")

# find base part

# over coordinates
        for coords in sorted(coords_dict):
    # write pos and source for each coordinate
            output.write(''.join(["pos:", str(coords), "(source)", "\n"]))
            name, sequence = fasta_dict[source].id, str(fasta_dict[source].seq)
            source_base = sequence[coords-1]
            output.write("genome: " + source + "\t" + "name: " + name + "\t" + "base: " + source_base + "(source)" + "\n")
    # list of difference for each coordinate
            dif = []
    # over destinations
            for dest in sorted(coords_dict[coords]):
                if dest == 'c':
                    name = "consensus"
                else:
                    name = fasta_dict[dest].id
        # find base when its on the - strand
                if coords_dict[coords][dest] < 0:
                    sequence = fasta_dict[dest].seq
                    base = Seq(sequence[((coords_dict[coords][dest])*-1)-1])
            # complement it to find base on the - strand
                    base = str(base.complement())
        # find base when its on the + strand
                else:
                    sequence = str(fasta_dict[dest].seq)
                    base = sequence[coords_dict[coords][dest]-1]
        # if base different: add in dif list
                if base != source_base:
                    dif.append(dest)
        # write output for one coordinate on one destination
                output.write("genome: " + dest + "\t" + "name: " + name + "\t" + "pos: " + str(coords_dict[coords][dest]) + "\t" + "base: " + base + "\n")
        # write dif list
            for item in dif:
                output.write(item + " ")
            output.write("\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--output_path", dest="output_p", help="path to output directory", required=True)
    parser.add_argument("-n", "--name", dest="output_name",
                        help="file prefix and sequence header for consensus FASTA / XFMA file", required=True)
    parser.add_argument("-c", "--consensus", dest="consensus_f", help="consensus FASTA file used in XMFA",
                        required=True)
    parser.add_argument("-i", "--index", dest="coord_f",
                        help="file with indices to map. First line: source_seq\tdest_seq[,dest_seq2,...] using \"c\" or sequence number. Then one coordinate per line. Coordinates are 1-based!")
    args = parser.parse_args()

    main()
