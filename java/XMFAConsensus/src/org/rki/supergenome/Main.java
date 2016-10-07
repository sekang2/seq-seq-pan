package org.rki.supergenome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Main {
		
	public static Map<String, List<String>> getParams(String[] args) {
		Map<String, List<String>> params = new HashMap<>();

		List<String> options = null;
		for (int i = 0; i < args.length; i++) {
			final String a = args[i];

			if (a.charAt(0) == '-') {
				if (a.length() < 2) {
					System.err.println("Error at argument " + a);
					return null;
				}

				options = new ArrayList<>();
				params.put(a.substring(1), options);
			}
			else if (options != null) {
				options.add(a);
			}
			else {
				System.err.println("Illegal parameter usage");
				return null;
			}
		}
		return params;
	}

	public static void main(String[] args) {
		Map<String, List<String>> params = getParams(args);
		
		XMFAReader reader = new XMFAReader();
		char[] spacer_c = new char[1000];
		for(int i=0; i<1000; i++)
			spacer_c[i] = 'N';
		String spacer = new String(spacer_c);
		
//		String infilename = "C:/temp/supergenome/add1_realign.xmfa";
//		String out_fasta = "C:/temp/supergenome/add1_consensus.fasta";
//		String out_fasta_idx = "C:/temp/supergenome/add1_consensus.fasta.idx";
//		String out_fasta_blocksep = "C:/temp/supergenome/add1_consensus.fasta.blockseparated.fasta";
//		String out_fasta_blocksep_idx = "C:/temp/supergenome/add1_consensus.fasta.blockseparated.idx";

		String infilename = "", out_fasta = "", out_fasta_idx = "", out_fasta_blocksep = "", out_fasta_blocksep_idx = "";
		
		try {
			infilename = params.get("-xmfa").get(0);
			out_fasta = params.get("-outfasta").get(0);
			out_fasta_idx = params.get("-outidx").get(0);
			out_fasta_blocksep = params.get("-outblocksepfasta").get(0);
			out_fasta_blocksep_idx = params.get("-outblocksepidx").get(0);
		} catch(NullPointerException e) {
			System.out.println("Please provide the following parameters:\n--xmfa (input xmfa file)\n--outfasta (output FASTA consensus)\n--outidx (idx file for output consensus with gap positions)\n--outblocksepfasta (output FASTA consensus with blocks separated by 1000 Ns)\n--outblocksepidx (output idx file for FASTA consensus with blocks separated by 1000 Ns)");
			System.exit(0);
		}

		String[] namepart_data = infilename.split("[\\/]");
		String namepart = namepart_data[namepart_data.length-1].split("_")[0];
		
		try {
			int numBlock = 0;
			int currentlen = 0;
			int blockseplen = 1;
			BufferedReader infile = new BufferedReader(new FileReader(new File(infilename)));
			BufferedWriter outfile_blocksep_fasta = new BufferedWriter(new FileWriter(new File(out_fasta_blocksep)));
			BufferedWriter outfile_blocksep_fasta_idx = new BufferedWriter(new FileWriter(new File(out_fasta_blocksep_idx)));
			BufferedWriter outfile_fasta = new BufferedWriter(new FileWriter(new File(out_fasta)));
			BufferedWriter outfile_fasta_idx = new BufferedWriter(new FileWriter(new File(out_fasta_idx)));
			reader.skipHead(infile);
			//TODO: Fix header
			outfile_blocksep_fasta.write(">" + namepart + ";0|" + infilename + "\n");
			outfile_fasta.write(">" + namepart + ";0|" + infilename + "\n");
			
			outfile_blocksep_fasta_idx.write("#Fasta\t" + out_fasta + "\n#XMFA\t" + infilename + "\n");
			outfile_fasta_idx.write("#Fasta\t" + out_fasta + "\n");
			outfile_fasta_idx.write("#XMFA\t" + infilename + "\n");
			outfile_fasta_idx.write(reader.getHeader() + "\n");
			while(infile.ready()) {
				numBlock += 1;
				LCB lcb = reader.readLCB(infile);
				String consensus = lcb.getConsensus();
//				System.out.print(lcb.numSequences);
//				System.out.print(": ");
//				System.out.println(consensus.length());
				if(currentlen != 0) {
					outfile_blocksep_fasta.write(spacer);
				}
				outfile_blocksep_fasta.write(consensus);
				outfile_fasta.write(consensus);
				if(currentlen != 0)
					outfile_blocksep_fasta_idx.write(";");
				outfile_blocksep_fasta_idx.write(new Integer(currentlen).toString());
				
				if(numBlock != 1)
					outfile_fasta_idx.write("\n");
				outfile_fasta_idx.write("b\t" + numBlock + "\t" + blockseplen + "\t");
				currentlen += consensus.length();
				blockseplen += consensus.length();
				currentlen += spacer.length();
				outfile_fasta_idx.write((blockseplen-1) + "\t+\n");
				outfile_fasta_idx.write(lcb.writeSeqBlocks());
			}
			infile.close();
			outfile_blocksep_fasta.close();
			outfile_blocksep_fasta_idx.close();
			outfile_fasta.close();
			outfile_fasta_idx.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
