/*
 * Quantification analysis based on Maxquant output
 */
import java.util.*;
import java.io.*;
public class Quant {
	public static void main(String[] args) {
	}
	
	//////////////////////////read raw output and store into an arraylist
	public static ArrayList<PepQuanDB> readFile(String filename) {
		ArrayList<PepQuanDB> apepQuanDB = null;
		Map<String, String> proteinMap = new HashMap<String, String>();
		PepQuanDB pepQuanDB;
		if (filename == null) return null;
		try {
			FileReader fr = new FileReader(filename);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			apepQuanDB = new ArrayList<PepQuanDB>();
			pepQuanDB = new PepQuanDB();
			while (line != null) {
				//System.out.println(line);
				String[] parseLine = line.split("\t");
				//System.out.println(parseLine.length);
				if (parseLine[0].equals("id")) {
					line = br.readLine();
					continue;
				}
				if (!proteinMap.containsKey(parseLine[11])){
					pepQuanDB = new PepQuanDB();
					apepQuanDB.add(pepQuanDB);
					proteinMap.put(parseLine[11], "TEST");
				}
				pepQuanDB.groupId = Integer.valueOf(parseLine[1]);
				pepQuanDB.proteinName = parseLine[11];
				PepDB pepDB = new PepDB();
				pepDB.pepSeq = parseLine[7];
				for (int i = 0; i < 9; i++)	pepDB.area[i] = Long.valueOf(parseLine[40 + i]);
				pepQuanDB.apepDB.add(pepDB);
				line = br.readLine();
			}
			System.out.println("The file reading is finished");
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		return apepQuanDB;
	}
	///////////////////////end of store raw output
	
	//calculate ratio for each peptide on per run basis; store all the measurement for one protein together
	public static ArrayList<ProDB> storePepRatio(ArrayList<PepQuanDB> aPepQuanDB) {
		ArrayList<ProDB> aProDB = new ArrayList<ProDB>();
		if (aPepQuanDB == null) return null;
		for(PepQuanDB pepQuanDB : aPepQuanDB) {
			ProDB proDB = new ProDB();
			proDB.proName = pepQuanDB.proteinName;
			//System.out.println(pepQuanDB.apepDB.size());
			//System.out.println(pepQuanDB.proteinName);
			//System.out.println(proDB.proName);
			for (PepDB pepDB : pepQuanDB.apepDB) {
				for (int i = 0; i < 3; i++) {//skip the measurement if a '0' appears as the area
					if (pepDB.area[i] != 0 && pepDB.area[i + 3] != 0) {
						double ratioML = MaxQuantMath.logE((double) pepDB.area[i + 3] / pepDB.area[i]);
						proDB.proRatioML.add(ratioML);
						//System.out.println(ratioML);
					}
					if (pepDB.area[i] != 0 && pepDB.area[i + 6] != 0) {
						double ratioUL = MaxQuantMath.logE((double) pepDB.area[i + 6] / pepDB.area[i]);
						proDB.proRatioUL.add(ratioUL);
						//System.out.println(ratioUL);
					}
				}
			}
			//System.out.println(proDB.proRatioML.size());
			aProDB.add(proDB);
		}
		System.out.println("Finished data storage");
		return aProDB;
	}
	///////////////////////////////end of store peptide ratios
	
	//////////////////////////////calculate the protein ratio and the stdev
	public static ArrayList<StatDB> storeProRatio(ArrayList<ProDB> aProDB) {
		ArrayList<StatDB> aStatDB = new ArrayList<StatDB>();
		if (aProDB == null) return null;
		for (ProDB proDB : aProDB) {
			StatDB statDB = new StatDB();
			statDB.proName = proDB.proName;
			//System.out.print(statDB.proName);
			//System.out.print(proDB.proRatioML.size() + "\t");
			//System.out.print(proDB.proRatioUL.size());
			//System.out.println();
			statDB.avgProRatioML = MaxQuantMath.avg(proDB.proRatioML);
			statDB.avgProRatioUL = MaxQuantMath.avg(proDB.proRatioUL);
			statDB.stdevProRatioML = MaxQuantMath.stdev(proDB.proRatioML);
			statDB.stdevProRatioUL = MaxQuantMath.stdev(proDB.proRatioUL);
			//statDB.avgProRatioMLlog = Math.log(statDB.avgProRatioML);
			//statDB.avgProRatioULlog = Math.log(statDB.avgProRatioUL);
			aStatDB.add(statDB);
		}
		System.out.println("The protien ratio calculation and storage is finished");
		return aStatDB;
	}
	////////////////////////////////end of protein ratio calculation and storage
	
	//////////////////////////// Print the result into a .txt file
	public static void printIntoFile(ArrayList<StatDB> aStatDB, String filename) {
		try {
			FileOutputStream fout = new FileOutputStream(filename);
			PrintStream ps = new PrintStream(fout);
			for (StatDB statDB : aStatDB) {
				String tmp = statDB.proName + "\t" + statDB.avgProRatioMLlog + "\t" + statDB.stdevProRatioML + "\t" + statDB.avgProRatioULlog + "\t" + statDB.stdevProRatioUL;
				ps.println(tmp);
			}
			System.out.println("The printing service is finished");
			ps.close();
			fout.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}
	/////////////////////////// end of file printing
	
	/*
	 * Read maxquant raw output, peptide view, all peptides.
	 * retain unique peptides, require each peptide has 2 out of 3 measurement in technical rep, 
	 * remove decoy and contamination hits, remove modified peptides
	 */
	public static void processRawOutput(String filename, String outputFilename) {
		if (filename == null) {
			System.out.println("No filename provided");
			return;
		}
		Map<String, PepQuanDB> proMap = new HashMap<String, PepQuanDB>();
		ArrayList<PepQuanDB> aPepQuanDB = new ArrayList<PepQuanDB>();
		try {
			FileReader fr = new FileReader(filename);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();			
			while (line != null) {
				//System.out.println(parseLine.length);
				//System.out.println(line);
				String[] parseLine = line.split("\t");
				//System.out.println(parseLine.length);
				if (parseLine.length > 50) {
					//System.out.println(line);
					line = br.readLine();
					continue;
				}
				if (!parseLine[13].equals("yes")) {
					line = br.readLine();
					continue;
				}
				if (!parseLine[6].equals("")) {
					line = br.readLine();
					continue;
				}
				//System.out.println(line);
				if (!proMap.containsKey(parseLine[11])) {
					PepQuanDB pepQuanDB = new PepQuanDB();
					pepQuanDB.proteinName = parseLine[11];
					pepQuanDB.groupId = Integer.valueOf(parseLine[1]);
					proMap.put(parseLine[11], pepQuanDB);
					aPepQuanDB.add(pepQuanDB);
				}
				PepQuanDB tmp = proMap.get(parseLine[11]);
				PepDB pepDB = new PepDB();
				pepDB.pepSeq = parseLine[7];
				for (int i = 0; i < 9; i++) {
					pepDB.area[i] = Long.valueOf(parseLine[41 + i]);
				}
				tmp.apepDB.add(pepDB);
				line = br.readLine();
			}
			System.out.println("The file reading is finished");
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		try {
			FileOutputStream fout = new FileOutputStream(outputFilename);
			PrintStream ps = new PrintStream(fout);
			for (PepQuanDB pepQuanDB : aPepQuanDB) {
				String printTmp = "";
				int pepCounter = 0;
				for (PepDB pepDB : pepQuanDB.apepDB) {
					int zeroCounter = 0;
					boolean lower = false;  //////////2 out of 3 runs in lower lobe have signal
					boolean middle = false; //////////2 out of 3 runs in middle lobe have signal
					boolean upper = false;  //////////2 out of 3 runs in the upper lobe have signal
					for (int i = 0; i < 3; i++) {
						if (pepDB.area[i] == 0) zeroCounter++;
					}
					if(zeroCounter <= 1)lower = true;
					zeroCounter = 0;
					for (int i = 3; i < 6; i++) {
						if (pepDB.area[i] == 0) zeroCounter++;
					}
					if (zeroCounter <= 1) middle = true;
					zeroCounter = 0;
					for (int i = 6; i < 9; i++) {
						if (pepDB.area[i] == 0) zeroCounter++;
					}
					if (zeroCounter <= 1) upper = true;
					if (lower && middle && upper) {
						pepCounter++;
						printTmp += pepQuanDB.groupId + "\t" + pepQuanDB.proteinName + "\t" + pepDB.pepSeq;
						for (int i = 0; i < 9; i++) printTmp += "\t" + pepDB.area[i];
						printTmp += "\n";
					}
				}
				if (pepCounter > 1) ps.print(printTmp);
			}
			System.out.println("The file writing is finished");
			ps.close();
			fout.close();
		} catch (IOException e){
			System.err.println(e.getMessage());
		}
		
	}
//////////////////end of the raw output analysis	
	
/*
 * read in the processed raw output files and do the following experiments,
 * take in the output from processRawOutput
 */
	public static void processWrapper(String filename, String outputFileName) {
		if (filename == "" || filename == null) return;
		Map<String, PepQuanDB> proMap = new HashMap<String, PepQuanDB>();
		ArrayList<PepQuanDB> aPepQuanDB = new ArrayList<PepQuanDB>();
		try {
			FileReader fr = new FileReader(filename);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				String[] parseLine = line.split("\t");
				if (parseLine[0].equals("id")) {
					line = br.readLine();
					continue;
				}
				if (!proMap.containsKey(parseLine[1])) {
					PepQuanDB pepQuanDB = new PepQuanDB();
					pepQuanDB.proteinName = parseLine[1];
					pepQuanDB.groupId = Integer.valueOf(parseLine[0]);
					proMap.put(parseLine[1], pepQuanDB);
					aPepQuanDB.add(pepQuanDB);
				}
				PepQuanDB buffer = proMap.get(parseLine[1]);
				PepDB bufferPepDB = new PepDB();
				bufferPepDB.pepSeq = parseLine[2];
				for (int i = 0; i < 9; i++) bufferPepDB.area[i] = Long.valueOf(parseLine[3 + i]);
				buffer.apepDB.add(bufferPepDB);
				line = br.readLine();
			}
			br.close();
			fr.close();
			System.out.println("The file reading is finished");
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		//ArrayList<ProDB> aProDB = storePepRatio(aPepQuanDB);
		ArrayList<StatDB> aStatDB = sumAreaforRatio(aPepQuanDB);
		//ArrayList<StatDB> aStatDB = storeProRatio(aProDB);
		//for (Double d : aProDB.get(0).proRatioML) System.out.println(d);
		printIntoFile(aStatDB, outputFileName);
	}
	/////////////////////////////////End of the wrapper
	
	//////////////////////////Another take ratio method, sum the area of all the peptides from one protein first, then take the ratio
	public static ArrayList<StatDB> sumAreaforRatio(ArrayList<PepQuanDB> aPepQuanDB) {
		if (aPepQuanDB == null) return null;
		ArrayList<StatDB> aStatDB = new ArrayList<StatDB>();
		for (PepQuanDB pepQuanDB : aPepQuanDB) {
			if (pepQuanDB.apepDB.size() < 2) continue;
			StatDB statDB = new StatDB();
			statDB.proName = pepQuanDB.proteinName;
			double[] areaTmp = new double[9];
			for (int i = 0; i < 9; i++) areaTmp[i] = 0;
			for (PepDB pepDB : pepQuanDB.apepDB) {
				for (int i = 0; i < 9; i++) {
					areaTmp[i] += pepDB.normalizedA[i];
				}
			}
			ArrayList<Double> ratioML = new ArrayList<Double>();
			ArrayList<Double> ratioUL = new ArrayList<Double>();
			for (int i = 0; i < 3;  i++) {
				if (areaTmp[i] != 0) {
					ratioML.add(MaxQuantMath.logE(areaTmp[i + 3] / areaTmp[i]));
					ratioUL.add(MaxQuantMath.logE(areaTmp[i + 6] / areaTmp[i]));
				}
			}
			statDB.avgProRatioMLlog = MaxQuantMath.avg(ratioML);
			statDB.stdevProRatioML = MaxQuantMath.stdev(ratioML);
			statDB.avgProRatioULlog = MaxQuantMath.avg(ratioUL);
			statDB.stdevProRatioUL = MaxQuantMath.stdev(ratioUL);
			aStatDB.add(statDB);
		}
		return aStatDB;
	}
	////////////////////////////////End of function, sum area first, then take the ratio
	
	///////////////////////////////Align pool samples and variance measurement within a pool to show the difference in the %stdev
	public static void comparePoolwithIsolates(String poolFilename, String isolatesFilename, String outputFilename) {
		if (poolFilename == null || poolFilename == "") {
			System.out.println("Pool filename missing");
			return;
		}
		if (isolatesFilename == null || isolatesFilename == "") {
			System.out.println("isolates filename missing");
			return;
		}
		Map<String, ProDB> proMap = new HashMap<String, ProDB>();
		ArrayList<ProDB> aProDB = new ArrayList<ProDB>();
		try {
			FileReader fr = new FileReader(poolFilename);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				String[] parseLine = line.split("\t");
				ProDB proDB = new ProDB();
				proDB.proName = parseLine[0].substring(0, 6);
				//System.out.println(proDB.proName);
				proDB.proRatioSumML = Double.valueOf(parseLine[1]);
				proDB.proRatioSumMLstdev = Math.abs(Double.valueOf(parseLine[2]) / proDB.proRatioSumML) * 100;
				proDB.proRatioSumUL = Double.valueOf(parseLine[3]);
				proDB.proRatioSumULstdev = Math.abs(Double.valueOf(parseLine[4]) / proDB.proRatioSumUL) * 100;
				//System.out.println(proDB.proRatioSumML + " " + proDB.proRatioSumMLstdev + " " + proDB.proRatioSumUL + " " + proDB.proRatioSumULstdev);
				proMap.put(proDB.proName, proDB);
				aProDB.add(proDB);
				line = br.readLine();
			}
			System.out.println("poolfile has been read");
			br.close();
			fr.close();
			FileReader fr1 = new FileReader(isolatesFilename);
			BufferedReader br1 = new BufferedReader(fr1);
			line = br1.readLine();
			while (line != null) {
				String[] parseLine = line.split("\t");
				//System.out.println(parseLine[0]);
				if (proMap.containsKey(parseLine[0])) {
					ProDB proDB = proMap.get(parseLine[0]);
					proDB.proRatioSum71 = Double.valueOf(parseLine[1]);
					proDB.proRatioSum71stdev = Math.abs(Double.valueOf(parseLine[2]) / proDB.proRatioSumML) * 100;
					proDB.proRatioSum81 = Double.valueOf(parseLine[3]);
					proDB.proRatioSum81stdev = Math.abs(Double.valueOf(parseLine[4]) / proDB.proRatioSumUL) * 100;
					//System.out.println(proDB.proRatioSum71 + " " + proDB.proRatioSum71stdev + " " + proDB.proRatioSum81 + " " + proDB.proRatioSum81stdev);
				}
				line = br1.readLine();
			}
			System.out.println("isolatefile has been read");
			br1.close();
			fr1.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		try {
			FileOutputStream fout = new FileOutputStream(outputFilename);
			PrintStream ps = new PrintStream(fout);
			//String printBuffer;
			for (ProDB proDB : aProDB) {
				if (proDB.proRatioSum71 * proDB.proRatioSum81 * proDB.proRatioSumML * proDB.proRatioSumUL != 0) {
					ps.println(proDB.proName + "\t" + 
							   proDB.proRatioSum71 + "\t" + 
							   proDB.proRatioSum71stdev + "%" + "\t" +
							   proDB.proRatioSum81 + "\t" +
							   proDB.proRatioSum81stdev + "%" + "\t" +
							   proDB.proRatioSumML + "\t" + 
							   proDB.proRatioSumMLstdev + "%" + "\t" +
							   proDB.proRatioSumUL + "\t" +
							   proDB.proRatioSumULstdev + "%" + "\t");
				}
			}
			ps.close();
			fout.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		System.out.println("The output file has been successfully created");
	}
///////////////////////////////end of the function

////////////////Store peptide information in a PepQuanDB
	public static ArrayList<PepDB> storePepQuanDB(String filename) {
		ArrayList<PepDB> aPepDB = new ArrayList<PepDB>();
		if (filename == null || filename == "") return null;
		try {
			FileReader fr = new FileReader(filename);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				String[] parseLine = line.split("\t");
				if (parseLine[0].equals("id")) {
					line = br.readLine();
					continue;
				}
				if (!parseLine[13].equals("yes")) {
					line = br.readLine();
					continue;
				}
				if (parseLine.length > 50) {
					line = br.readLine();
					continue;
				}
				PepDB pepDB = new PepDB();
				pepDB.pepSeq = parseLine[7];
				if(parseLine[11].substring(0, 1).equals("\"")) pepDB.proName = parseLine[11].substring(1, 7);
				else pepDB.proName = parseLine[11].substring(0, 6);
				for (int i = 0; i < 9;  i++) {
					pepDB.area[i] = Long.valueOf(parseLine[i + 41]);
				}
				line = br.readLine();
				aPepDB.add(pepDB);
			}
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		System.out.println(filename + " has been successfully read");
		return aPepDB;
	}
///////////////End of the function

/////////////// filter the peptides by requiring 9 out of 9 presence and lower than 50% RSD within technical replicates
	public static ArrayList<PepDB> pepFilter(ArrayList<PepDB> aPepDB, String outputfilename, double percent) {
		if (aPepDB == null) return null;
		ArrayList<PepDB> bPepDB = new ArrayList<PepDB>();
		//System.out.println(aPepDB.size());
		/*
		 * Calculate the normalization ratio for each sample
		 */
		long[] sum = new long[9];
		for (PepDB pepDB : aPepDB) {
			for (int i = 0; i < 9; i++) {
				sum[i] += pepDB.area[i];
			}
		}
		double[] ratio = new double[9];
		for (int i = 0; i < 9; i++) {
			ratio[i] = (double) sum[0] / sum[i];
			//ratio[i] = 1;
		}
		////////////// end of the normalization ratio calculation
		for (PepDB pepDB : aPepDB) {
			boolean allPresence = true;
			boolean isRSDsmall = true;
			//System.out.println(pepDB.pepSeq);
			for (int i = 0; i < 9; i++) {
				if (pepDB.area[i] == 0) {
					allPresence = false;
					break;
				}
			}
			if (allPresence) {
				for (int i = 0; i < 9;) {
					double avg = (pepDB.area[i] * ratio[i] + pepDB.area[i + 1] * ratio[i + 1] + pepDB.area[i + 2] * ratio[i + 2]) / 3;
					double stdev = (double)(pepDB.area[i] * pepDB.area[i] * ratio[i] * ratio[i]
											+ pepDB.area[i + 1] * pepDB.area[i + 1] * ratio[i + 1] * ratio[i + 1] 
											+ pepDB.area[i + 2] * pepDB.area[i + 2] * ratio[i + 2] * ratio[i + 2]) / 3 - avg * avg;
					//System.out.println(avg + " " + Math.sqrt(stdev));
					//System.out.println(Math.sqrt(Math.abs(stdev / avg)));
					if (Math.sqrt(stdev) / avg > percent) {
						isRSDsmall = false;
						break;
					}
					i += 3;
				}
			}
			if (isRSDsmall && allPresence) bPepDB.add(pepDB); 
		}
		for (PepDB pepDB : bPepDB) {
			for (int i = 0; i < 9; i++) {
				pepDB.normalizedA[i] = pepDB.area[i] * ratio[i];
			}
		}
		/*for (PepDB pepDB : bPepDB) {
			for (int i = 0; i < 9; i++) {
				sum[i] += pepDB.area[i];
			}
		}
		for (int i = 0; i < 9; i++) {
			ratio[i] = (double) sum[0] / sum[i];
			//ratio[i] = 1;
		}
		ArrayList<PepDB> cPepDB = new ArrayList<PepDB>();
		for (PepDB pepDB : bPepDB) {
			boolean isRSDsmall = true;
			for (int i = 0; i < 9;) {
				double avg = (pepDB.area[i] * ratio[i] + pepDB.area[i + 1] * ratio[i + 1] + pepDB.area[i + 2] * ratio[i + 2]) / 3;
				double stdev = (double)(pepDB.area[i] * pepDB.area[i] * ratio[i] * ratio[i]
										+ pepDB.area[i + 1] * pepDB.area[i + 1] * ratio[i + 1] * ratio[i + 1] 
										+ pepDB.area[i + 2] * pepDB.area[i + 2] * ratio[i + 2] * ratio[i + 2]) / 3 - avg * avg;
				//System.out.println(avg + " " + Math.sqrt(stdev));
				//System.out.println(Math.sqrt(Math.abs(stdev / avg)));
				if (Math.sqrt(stdev) / avg > 0.5) {
					isRSDsmall = false;
					break;
				}
				i += 3;
			}
			if (isRSDsmall) cPepDB.add(pepDB);
		}*/
		if (outputfilename != null && outputfilename != "") {
			try {
				FileOutputStream fout = new FileOutputStream(outputfilename);
				PrintStream ps = new PrintStream(fout);
				for (PepDB pepDB : bPepDB) {
				//for (PepDB pepDB : cPepDB) {
					ps.print(pepDB.pepSeq + "\t" +
							 pepDB.proName + "\t");
					for (int i = 0; i < 9; i++) {
						ps.print(pepDB.normalizedA[i] + "\t");
					}
					ps.println();
				}
				ps.close();
				fout.close();
				System.out.println(outputfilename + " has been successfuly created!");
			} catch (IOException e) {
				System.out.println(e.getMessage());
			}
		}		
		return bPepDB;
		//return cPepDB;
	}
/////////////// End of the function
	
///////////////store the apepDb into a pepQuanDB array list
	public static ArrayList<PepQuanDB> storePepDBintoPepQuanDB(ArrayList<PepDB> aPepDB) {
		 ArrayList<PepQuanDB> aPepQuanDB = new ArrayList<PepQuanDB>();
		 Map<String, PepQuanDB> proMap = new HashMap<String, PepQuanDB>();
		 for (PepDB pepDB : aPepDB) {
			 if (!proMap.containsKey(pepDB.proName)) {
				 PepQuanDB pepQuanDB = new PepQuanDB();
				 pepQuanDB.proteinName = pepDB.proName;
				 pepQuanDB.apepDB.add(pepDB);
				 proMap.put(pepDB.proName, pepQuanDB);
				 aPepQuanDB.add(pepQuanDB);
			 } else {
				 PepQuanDB tmp = proMap.get(pepDB.proName);
				 tmp.apepDB.add(pepDB);
			 }
		 }
		 return aPepQuanDB;
	}
////////////////////finish implementation of the function
	
///////////////////function to return all the pep ratios from one set of technical replicates, write that results into a file
	public static void returnAllRatiosfromTechRep(ArrayList<PepDB> aPepDB, String outputFilename) {
		if (aPepDB == null) {
			System.out.println("Input database is null");
			return;
		}
		if (outputFilename == "" || outputFilename == null) {
			System.out.println("Please give a file name for the output");
		}
		try {
			FileOutputStream fout = new FileOutputStream(outputFilename);
			PrintStream ps = new PrintStream(fout);
			for (PepDB pepDB : aPepDB) {
				ps.println(pepDB.pepSeq + "\t" + 
						   pepDB.proName + "\t" +
						   pepDB.normalizedA[1] / pepDB.normalizedA[0] + "\t" +
						   pepDB.normalizedA[2] / pepDB.normalizedA[0] + "\t" +
						   pepDB.normalizedA[3] / pepDB.normalizedA[0] + "\t" +
						   pepDB.normalizedA[4] / pepDB.normalizedA[0] + "\t" +
						   pepDB.normalizedA[5] / pepDB.normalizedA[0] + "\t" +
						   pepDB.normalizedA[6] / pepDB.normalizedA[0] + "\t" +
						   pepDB.normalizedA[7] / pepDB.normalizedA[0] + "\t" +
						   pepDB.normalizedA[8] / pepDB.normalizedA[0] + "\t");
			}
			ps.close();
			fout.close();
			System.out.println(outputFilename + " has been successfully created!");
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}
///////////////////end of the function

//////////// read in peptide ratios, convert peptide ratio into protein ratio and output them
	public static void pepDBToProRatioperRun(ArrayList<PepQuanDB> aPepQuanDB, String outputFilename) {
		if (aPepQuanDB == null) {
			System.out.println("Please input the database");
			return;
		}
		if (outputFilename == null || outputFilename == "") {
			System.out.println("Please give an output file name");
			return;
		}
		try {
			FileOutputStream fout = new FileOutputStream(outputFilename);
			PrintStream ps = new PrintStream(fout);
			for (PepQuanDB pepQuanDB : aPepQuanDB) {
				long[] areaTmp = new long[9];
				if (pepQuanDB.apepDB.size() < 2) continue;
				for (PepDB pepDB : pepQuanDB.apepDB) {
					for (int i = 0; i < 9; i++) areaTmp[i] += pepDB.normalizedA[i];
				}
				String buffer = "";
				for (int i = 0; i < 9; i++) buffer = buffer + Double.toString(Math.log((double) 3 * areaTmp[i] / (areaTmp[0] + areaTmp[1] + areaTmp[2]))) + "\t";
				ps.println(pepQuanDB.proteinName + "\t" + buffer);
			}
			ps.close();
			fout.close();
			System.out.println("The file writing job is finished");
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		
	}
//////////// end of the function
}
