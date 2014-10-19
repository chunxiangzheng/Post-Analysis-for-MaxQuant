/*
 * file parsing and related processing methods for the heat map analysis
 */
import java.util.*;
import java.io.*;
public class HeatMapProcess {
	public static void main(String[] args) {
		fileReadWrite("HeatMapInput.txt", "HeatMapReady.txt");
	}
	/*
	 * find the proteins which are common in all the datasets and output them in a right format cluster can take
	 */
	public static void fileReadWrite(String input, String output) {
		if (input == null || input == "") {
			System.out.println("Input file name is needed");
			return;
		}
		if (output == null || output == "") {
			System.out.println("Output file name is needed");
			return;
		}
		Map<String, ProRatio> proMap = new HashMap<String, ProRatio>();
		ArrayList<ProRatio> aProRatio = new ArrayList<ProRatio>();
		try {
			FileReader fr = new FileReader(input);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				String[] parseLine = line.split("\t");
				if (parseLine[0].equals("Pro")) {
					line = br.readLine();
					continue;
				}
				if (parseLine.length == 7) {
					if (!proMap.containsKey(parseLine[4])) {
						ProRatio proRatio1 = new ProRatio();
						proRatio1.proname = parseLine[4];
						aProRatio.add(proRatio1);
						proMap.put(parseLine[4], proRatio1);
					}
					ProRatio proRatio1 = proMap.get(parseLine[4]);
					proRatio1.ratio[3] = parseLine[5];
					proRatio1.ratio[4] = parseLine[6];
					
				}
				if (!proMap.containsKey(parseLine[0])){
					ProRatio proRatio = new ProRatio();
					proRatio.proname = parseLine[0];
					aProRatio.add(proRatio);
					proMap.put(parseLine[0], proRatio);
				}
				ProRatio proRatio = proMap.get(parseLine[0]);
				for (int i = 1; i < 4; i++) {
					proRatio.ratio[i - 1] = parseLine[i];
				}
				line = br.readLine();
			}
			System.out.println("The file have been succcesfully read in.");
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		try {
			FileOutputStream fout = new FileOutputStream(output);
			PrintStream ps = new PrintStream(fout);
			for (ProRatio proRatio : aProRatio) {
				boolean isAll = true;
				for (int i = 0; i < 5; i++) {
					if (proRatio.ratio[i] == null) isAll = false;
				}
				if (isAll) {
					ps.println(proRatio.proname + "\t" +
							   proRatio.ratio[0] + "\t" +
							   proRatio.ratio[1] + "\t" +
							   proRatio.ratio[2] + "\t" +
							   proRatio.ratio[3] + "\t" +
							   proRatio.ratio[4]);
				}
			}
			System.out.println("The printing service is finished");
			ps.close();
			fout.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}
}
class ProRatio{
	String proname;
	String[] ratio;
	public ProRatio() {
		ratio = new String[5];
	}
}
