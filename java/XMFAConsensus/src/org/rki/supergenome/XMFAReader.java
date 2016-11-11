package org.rki.supergenome;

import java.io.BufferedReader;
import java.io.IOException;

public class XMFAReader {
	StringBuilder header = new StringBuilder();
	
	public XMFAReader() {
	}
	
	public void skipHead(BufferedReader reader) throws IOException {
		String line;
		reader.mark(2048);
		while((line=reader.readLine()) != null) {
			char head = line.charAt(0);
			if(line.startsWith("#Sequence")) {
				header.append(line);
				header.append("\n");
			}
			if(head != '#') {
				reader.reset();
				return;
			}
			reader.mark(2048);
		}		
	}
	
	public String getHeader() {
		return header.toString();
	}
	
	public LCB readLCB(BufferedReader reader) throws IOException {
		LCB res = new LCB();
		String line;
		while((line=reader.readLine()) != null) {
			if(!line.isEmpty()){
				char head = line.charAt(0);
				if(head == '=') {
					res.makeConsensus();
					break;
				}
				else if(head == '>')
					res.startNewSequence(line);
				else
					res.addLine(line.trim());
			}
		}
		return res;
	}
}
