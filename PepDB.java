class PepDB {
	public String pepSeq; //peptide sequence
	public String proName; // protein name
	public long[] area; // peptide area, 0-2 lower lobe 3-5 middle lobe 6-8 upper lobe
	public double[] normalizedA;
	public PepDB(){
		area = new long[9];
		normalizedA = new double[9];
	}
}
