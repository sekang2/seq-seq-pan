package org.rki.supergenome;

import java.util.LinkedList;

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
	public static char[] REVTABLE = {'-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

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
	private char[] consensus = null;
	private String consensusString = null;
	private Sequence currentSequence = null;
	
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
		for(int i=0; i<consensus.length; i++) {
			consensus[i] = REVTABLE[consensus[i]];
		}
		consensusString = new String(consensus);
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
		if(length == 0 && currentSequenceString != null) {
			consensus = currentSequenceString.toString().toCharArray();
			length = consensus.length;
			GapStretch gap = null;
			int lastGap = -2;
			for(int i=0; i<length; i++) {
				if(consensus[i] == '-') {
					if(gap == null) {
						gap = new GapStretch(i);
					}
					else if(lastGap != i-1) {
						gap.to = lastGap+1;
						sequences.getLast().addGap(gap);
						gap = new GapStretch(i);
					}
					lastGap = i;
				}
				else if(consensus[i] < TABLEOFFSET || consensus[i] > BASETABLE.length + TABLEOFFSET) {
					continue;
				}
				consensus[i] = BASETABLE[consensus[i] - TABLEOFFSET];
			}
			if(gap != null) {
				gap.to = lastGap+1;
				sequences.getLast().addGap(gap);				
			}
		}
		else if(currentSequenceString != null && currentSequenceString.length() != 0) {
			char[] seq = currentSequenceString.toString().toCharArray();
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
				}
				else if(seq[i] < TABLEOFFSET || seq[i] > BASETABLE.length + TABLEOFFSET) {
					continue;
				}
				consensus[i] |= BASETABLE[seq[i] - TABLEOFFSET];
			}
			if(gap != null) {
				gap.to = lastGap+1;
				sequences.getLast().addGap(gap);				
			}
		}
	}
	
}
