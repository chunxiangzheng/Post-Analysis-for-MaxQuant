/*
 * store protein name, each peptide ratio for each protein
 */
import java.util.*;
public class ProDB {
	String proName;
	double proRatioSumML;
	double proRatioSumUL;
	double proRatioSum71;
	double proRatioSum81;
	double proRatioSumMLstdev;
	double proRatioSumULstdev;
	double proRatioSum71stdev;
	double proRatioSum81stdev;
	ArrayList<Double> proRatioML;
	ArrayList<Double> proRatioUL;
	public ProDB(){
		proRatioML = new ArrayList<Double>();
		proRatioUL = new ArrayList<Double>();
	}
}
