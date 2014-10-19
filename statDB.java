/*
 * calculate and store protein ratio average and the stdev of the ratios
 */
public class StatDB {
	String proName; // protein name
	double avgProRatioML; // average protein ratio of all the peptide ratio for this protein (middle/lower lobe)
	double avgProRatioUL; // average protein ratio of all the peptide ratio for this protein (upper/lower lobe)
	double stdevProRatioML; // std deviation for the ratios ML
	double stdevProRatioUL; //std deviation for the ratios UL
	double avgProRatioMLlog;
	double avgProRatioULlog;
	public StatDB(){
	}
}
