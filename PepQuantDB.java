/*
 * quantification database generated from maxquant raw data
 */
import java.util.*;
public class PepQuanDB {
	public int groupId; //groupID from maxquant raw output
	public String proteinName; // protein name
	public ArrayList<PepDB> apepDB; //store peptide area information
	public long[] proArea;
	public PepQuanDB() {
		apepDB = new ArrayList<PepDB>();
		proArea = new long[9];
	}
}
