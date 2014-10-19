/*
 * store all the mathematical operations and functions here
 */
import java.util.*;
public class MaxQuantMath {
	public static double stdev(ArrayList<Double> al) {
		double sqrSum = 0;
		int count = al.size();
		if (count == 0) return 0;
		double avg = avg(al);
		for (Double a : al) sqrSum += a * a; 
		sqrSum = sqrSum/count - avg * avg;
		return  Math.sqrt(sqrSum); 
	}
	public static double avg(ArrayList<Double> al){
		double sum = 0;
		int count;
		count = al.size();
		if (count == 0) return 0;
		for (double a : al) {
			sum += a;
		}
		return sum/count;
	}
	public static double logE(double d) {
		return Math.log(d);
	}
}
