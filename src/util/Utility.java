package util;

import java.util.ArrayList;
import java.util.List;

public class Utility {

	
	public static double[] calcHistogram(double[] data, double min, double max, int numBins) {
		  final double[] result = new double[numBins];
		  final double binSize = (max - min)/numBins;

		  for (double d : data) {
		    int bin = (int) ((d - min) / binSize);
		    if (bin < 0) { /* this data is smaller than min */ }
		    else if (bin >= numBins) { /* this data point is bigger than max */ }
		    else {
		      result[bin] += 1.0;
		    }
		  }
		  /*double sum = 0;
		  for(double element: result)
			  sum += element;
			  
		  for(int i=0; i<result.length; i++){
			  result[i] = result[i] / ((max-min)/numBins*sum);
		  }
		  */
		  //hist = hist / ((((xmax-xmin)/no_bins)*hist.sum()))
		  
		  
		  return result;
		}
	
	public static double[] linearSpread(double min, double max, int numBins) {
		  final double[] result = new double[numBins];
		  final double step = (max-min)/ numBins;
		  
		  for(int i=0; i<numBins; i++){
			  result[i] = min + i*step;
		  }
		  return result;
		}
	
	// Order the result by id
	public static double[] orderResult(List<Double[]> resultRaw) {
		
		int index;
		double[] resultOuput = new double[resultRaw.size()];
		
		//System.out.println(resultRaw.size());
		
		for(int i=0 ; i<resultRaw.size();i++){
			index = (int) resultRaw.get(i)[0].intValue();
			resultOuput[index] = resultRaw.get(i)[1];
		}
		return resultOuput;
	}
	
	//
	// Average over regular and inverted sample pairs
	public static double[] processResult(List<Double[]> resultRaw) {
		
		double[] resultOrdered = orderResult(resultRaw);
		double[] resultOuput = new double[resultRaw.size()/2];
		
		for(int i=0 ; i<resultRaw.size()/2;i++){
			resultOuput[i] = (resultOrdered[2*i] + resultOrdered[2*i+1])/2;
		}
		return resultOuput;
	}
	
}
