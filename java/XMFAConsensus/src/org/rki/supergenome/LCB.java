package org.rki.supergenome;

import java.util.LinkedList;
import java.util.Random;

public class LCB {
	public static char BASE_A = 0b0001;
	public static char BASE_C = 0b0010;
	public static char BASE_G = 0b0100;
	public static char BASE_T = 0b1000;
	public static char BASE_R = (char)(BASE_A | BASE_G);
	public static char BASE_Y = (char)(BASE_C | BASE_T);
	public static char BASE_S = (char)(BASE_G | BASE_C);
	public static char BASE_W = (char)(BASE_A | BASE_T);
	public static char BASE_K = (char)(BASE_G | BASE_T);
	public static char BASE_M = (char)(BASE_A | BASE_C);
	public static char BASE_B = (char)(BASE_C | BASE_G | BASE_T);
	public static char BASE_D = (char)(BASE_A | BASE_G | BASE_T);
	public static char BASE_H = (char)(BASE_A | BASE_C | BASE_T);
	public static char BASE_V = (char)(BASE_A | BASE_C | BASE_G);
	public static char BASE_N = 0b1111;
	public static char BASE_GAP = 0b0000;
	public static int TABLEOFFSET = 45;
	public static char[] BASETABLE = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //[0]: '-', -> offset: 45
			BASE_A, BASE_B, BASE_C, BASE_D, 0, 0, BASE_G, BASE_H, 0, 0, BASE_K, 0, BASE_M, BASE_N, 0, 0, 0, BASE_R, BASE_S, BASE_T, 0, BASE_V, BASE_W, 0, BASE_Y, 0, //[20]: 'A' 
			0, 0, 0, 0, 0, 0, //[46]-[51]: '[' - '`'
			BASE_A, BASE_B, BASE_C, BASE_D, 0, 0, BASE_G, BASE_H, 0, 0, BASE_K, 0, BASE_M, BASE_N, 0, 0, 0, BASE_R, BASE_S, BASE_T, 0, BASE_V, BASE_W, 0, BASE_Y, 0}; //[52]: 'a'
//	public static char[] REVTABLE = {'-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

	
	public static int[] FACTOR = {0, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0};
	public static char[] UNAMB = {BASE_GAP, BASE_A, BASE_C, BASE_G, BASE_T};
	public static int UNAMB_LEN = 5;
	
	public static char[] REVTABLE = {'-', 'A', 'C', 'G', 'T'};
	
	private class GapStretch {
		public int from;
		public int to;
		
		public GapStretch(int from) {
			this.from = from;
			this.to = from;
		}
		
		public String toString() {
			return from + "-" + to;
		}
	}
	
	private class Sequence {
		private String num;
		private String start;
		private String end;
		private String strand;
		private LinkedList<GapStretch> gaps = new LinkedList<GapStretch>();
		
		public Sequence(String num, String start, String end, String strand) {
			this.num = num;
			this.start = start;
			this.end = end;
			this.strand = strand;
		}
		
		public void addGap(GapStretch gap) {
			gaps.add(gap);
		}
		
		public String toString() {
			StringBuilder res = new StringBuilder();
			res.append("s\t");
			res.append(num);
			res.append("\t");
			res.append(start);
			res.append("\t");
			res.append(end);
			res.append("\t");
			res.append(strand);
			res.append("\t");
			boolean firstGap = true;
			for(GapStretch gap : gaps) {
				if(!firstGap)
					res.append(";");
				res.append(gap);
				firstGap=false;
			}
			return res.toString();
		}
		
	}
	
	public int numSequences = 0;

	private StringBuilder currentSequenceString = null;
	private int length = 0;
	private LinkedList<Sequence> sequences = new LinkedList<Sequence>();
	private int[][] consensus = null;
	private String consensusString = null;
	private Sequence currentSequence = null;
	private Random random = new Random(23121989);
	
	public void startNewSequence(String defline) {
		if(currentSequence != null)
			sequences.add(currentSequence);
		String[] data = defline.split(" ");
		currentSequence = new Sequence(data[1].split(":")[0], data[1].split(":")[1].split("-")[0], data[1].split(":")[1].split("-")[1], data[2]);
		numSequences ++;
		if(currentSequenceString != null)
			addSequenceToConsensus();
		currentSequenceString = new StringBuilder();
	}
	
	public void addLine(String line) {
		currentSequenceString.append(line);
	}

	public void makeConsensus() {
		sequences.add(currentSequence);
		addSequenceToConsensus();
		char[] seq = new char[length];
		for(int i=0; i<length; i++) {
			/* loop through unambiguous base array and find max:
			 * if max is 0 use "N"
			 * if max is greater than 0 and current value is equal, change index based on random number
			 * if all values are the same use "N" 
			 * if any base is present, do not use gap 
			 */
			int max_idx = 0;
			int max_val = 0;
			int equal = 0;
			for(int j=0; j < UNAMB_LEN; j++){
				int cur_val = consensus[i][j];
				
				// current value is bigger or current max is a gap
				if(cur_val > max_val | (cur_val > 0 & max_idx == 0 )){ 
					max_val = cur_val;
					max_idx = j;
					equal = 0;
				}else if(cur_val == max_val){
					
					// change base based on random number
					if(random.nextBoolean()){
						max_idx = j;
					}
					
					equal += 1;
				}
			}
			
			// no max value or all values of bases the same
			if(max_val == 0 || equal == (UNAMB_LEN - 2)){
				seq[i] = 'N';
			}else{
				seq[i] = REVTABLE[max_idx];
			}
		}
		consensusString = new String(seq);
	}
	
	public String getConsensus() {
		if(consensusString == null)
			makeConsensus();
		return consensusString;
	}

	public String writeSeqBlocks() {
		StringBuilder res = new StringBuilder();
		for(Sequence seq : sequences) {
			res.append(seq);
			res.append("\n");
		}
		return res.toString();
	}
	
	private void addSequenceToConsensus() {
		if(currentSequenceString != null && currentSequenceString.length() != 0){
			
			char[] seq = currentSequenceString.toString().toCharArray();
		
			if(length == 0 ) {
				length = seq.length;
				consensus = new int[length][UNAMB_LEN];
			}
			
			
			GapStretch gap = null;
			int lastGap = -2;
			for(int i=0; i<length; i++) {
				if(seq[i] == '-') {
					if(gap == null) {
						gap = new GapStretch(i);
					}
					else if(lastGap != i-1) {
						gap.to = lastGap+1;
						sequences.getLast().addGap(gap);
						gap = new GapStretch(i);
					}
					lastGap = i;
					consensus[i][0] += 1;
				}
				else if(seq[i] < TABLEOFFSET || seq[i] > BASETABLE.length + TABLEOFFSET) {
					continue;
				}else{
					char bin_rep = BASETABLE[seq[i] - TABLEOFFSET];
					
					// loop through unambigious bases
					// add factor (3, 2, 1 ) to positions in consensus array
					for (int c=1; c < UNAMB_LEN; c++){
						if((bin_rep | UNAMB[c]) == bin_rep){
							consensus[i][c] += FACTOR[bin_rep];
						}
					}
				}
				
			}
			if(gap != null) {
				gap.to = lastGap+1;
				sequences.getLast().addGap(gap);				
			}
		}
	}
	
}
