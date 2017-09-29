import argparse
import os

from seqseqpan.io import Parser
from seqseqpan.mapper import Mapper
from Bio import SeqIO


def main():
    global args
    parser = Parser()
    mapper = Mapper()

    sparse_align, sparse_consensus = parser.parse_consensus_index(args.reference_f)
    ref_fasta = SeqIO.read(open(args.reference_f), 'fasta')
    alt_fasta = SeqIO.read(open(sparse_align.genomes[int(args.alternative_f)].file_path), 'fasta')
    source = "c"
    destination = args.alternative_f
    coords_ref = list(range(1, len(ref_fasta.seq)+1))
    # coords_ref = list(range(1,100000))
    coords_dict = mapper.map_coordinates(sparse_align, sparse_consensus, source, destination, coords_ref)

    with open(os.path.abspath(args.output_f + ".vcf"), "w") as output:
        output.write("##fileformat=VCFv4.2" + "\n")
        output.write("##Filter=<ID=150,Description=\"deletion longer then 150\">" + "\n")
        output.write("##INFO=<ID=TY,Number=1,Type=String,Description=\"Type\">" + "\n")
        output.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" + "\n")
        output.write("##reference=" + args.reference_f + "\n")
        output.write("#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\t" + alt_fasta.name + "\n")

        counter = 0
        del_length = 0
        ref_bases = ""
        for coords in sorted(coords_dict):

            # deletion
            if coords_dict[coords] == {}:
                del_length += 1
                ref_bases += ref_fasta.seq[coords-1]
                # end of ref or end of del
                if coords == len(ref_fasta.seq) or coords_dict[coords+1] != {}:
                    # del at pos 1
                    if coords-del_length == 0:
                        if del_length > 150:
                            output.write(ref_fasta.name + "\t" + str(1) + "\t" + "." + "\t" + ref_bases + ref_fasta.seq[coords] + "\t" + alt_fasta.seq[coords_dict[coords+1][args.alternative_f]-1] + "\t" + "." + "\t" + "150" + "\t" + "TY=DEL" + "\t" + "GT" + "\t" + "1" + "\n")
                        else:
                            output.write(ref_fasta.name + "\t" + str(1) + "\t" + "." + "\t" + ref_bases + ref_fasta.seq[coords] + "\t" + alt_fasta.seq[coords_dict[coords+1][args.alternative_f]-1] + "\t" + "." + "\t" + "PASS" + "\t" + "TY=DEL" + "\t" + "GT" + "\t" + "1" + "\n")
                    else:
                        if del_length > 150:
                            output.write(ref_fasta.name + "\t" + str(coords-del_length) + "\t" + "." + "\t" + ref_fasta.seq[coords-del_length-1] + ref_bases + "\t" + alt_fasta.seq[coords_dict[coords-del_length][args.alternative_f]-1] + "\t" + "." + "\t" + "150" + "\t" + "TY=DEL" + "\t" + "GT" + "\t" + "1" + "\n")
                        else:
                            output.write(ref_fasta.name + "\t" + str(coords-del_length) + "\t" + "." + "\t" + ref_fasta.seq[coords-del_length-1] + ref_bases + "\t" + alt_fasta.seq[coords_dict[coords-del_length][args.alternative_f]-1] + "\t" + "." + "\t" + "PASS" + "\t" + "TY=DEL" + "\t" + "GT" + "\t" + "1" + "\n")
                    counter += del_length
                    del_length = 0
                    ref_bases = ""

            # match
            elif ref_fasta.seq[coords-1] == alt_fasta.seq[coords_dict[coords][args.alternative_f]-1]:
                counter += 1

            # SNP
            elif ref_fasta.seq[coords-1] != alt_fasta.seq[coords_dict[coords][args.alternative_f]-1]:
                output.write(ref_fasta.name + "\t" + str(coords) + "\t" + "." + "\t" + ref_fasta.seq[coords-1] + "\t" + alt_fasta.seq[coords_dict[coords][args.alternative_f]-1] + "\t" + "." + "\t" + "PASS" + "\t" + "TY=SNP" + "\t" + "GT" + "\t" + "1" + "\n")
                counter += 1

            else:
                print("Fehler1:", coords)
        if counter == len(coords_ref):
            print("check")
        else:
            print("Fehler2")


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-r", "--reference", dest="reference_f", help="pangenome FASTA file as reference", required=True)
    arg_parser.add_argument("-a", "--alternative", dest="alternative_f", help="number of alternative genome", required=True)
    arg_parser.add_argument("-o", "--output", dest="output_f", help="VCF output file", required=True)
    args = arg_parser.parse_args()

    main()
