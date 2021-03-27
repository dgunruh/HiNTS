package routines;
import java.io.File;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
//import com.sun.corba.se.spi.orbutil.fsm.Guard.Result;

import classes.Nanoparticle;
import classes.Sample;
import util.Configuration;
import util.Constants;
//import org.apache.commons.math3.*;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

public class Setup {
	
	// Don't really need this method...
	public static Sample setupSample(Map<String,Object> params) throws Exception{
		
		Sample sample = new Sample(params);
		
		return sample;
	}
	
	// Should be called in the constructor for Sample
	public static Nanoparticle[] setupNanoparticles(Sample sample)  {
		
		List<String> line;
		double xcoord, ycoord, zcoord, diameter;
		int id;
		
		Nanoparticle[] nanoparticles = new Nanoparticle[sample.nnanops];
		List<String> lines;
		try {
			lines = Files.readLines(new File(sample.getFilename()), Charsets.UTF_8);
			for(int i=0; i<sample.nnanops; i++){
				//System.out.println(lines);
				//line = Arrays.asList(lines.get(i+4).split(" ")); Moule style
				if(!sample.necking) {
					line = Arrays.asList(lines.get(i+7).split(" "));
				}
				else {
					line = Arrays.asList(lines.get(i+7).split(" "));
				}
				
				//id = Integer.valueOf(line.get(0));
				/*xcoord = Double.valueOf(line.get(3));
				ycoord = Double.valueOf(line.get(4));
				zcoord = Double.valueOf(line.get(5));*/
				
				xcoord = Double.valueOf(line.get(3)); //5 for Moule
				ycoord = Double.valueOf(line.get(4)); //3 for Moule
				zcoord = Double.valueOf(line.get(5)); //4 for Moule
				
				if(sample.necking){
					xcoord = Double.valueOf(line.get(3));
					ycoord = Double.valueOf(line.get(4));
					zcoord = Double.valueOf(line.get(5));
				}
				
				if (sample.latticeStructure == "Moule"){
					xcoord = Double.valueOf(line.get(5)); //5 for Moule
					ycoord = Double.valueOf(line.get(3)); //3 for Moule
					zcoord = Double.valueOf(line.get(4)); //4 for Moule
				}
				diameter = Double.valueOf(line.get(6)); // in nm
				if(sample.necking) {
					diameter = Double.valueOf(line.get(6));
				}
				//System.out.println("diameter is: " + diameter);
				//diameter = 5.0;
				nanoparticles[i] = new Nanoparticle(xcoord, ycoord, zcoord, diameter, i, sample);
				sample.nanoparticle_idList.put(i, i);
			}
			return nanoparticles;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
		return nanoparticles;

	}
	
	
	/**
	@param nanoparticles and sample
	@return return the FWHM in ry unit
	@throws what kind of exception does this method throw
	*/
	public static double getFWHM(Nanoparticle[] nanoparticles, Sample sample){
		double FWHM;
		
		//Start by creating a histogram of the energy levels
		int bins = 50;
		double[] levels_eV = new double[sample.nnanops];
		double[] levels = new double[sample.nnanops];
		for(int i=0; i<sample.nnanops; i++){
			levels[i] = nanoparticles[i].getCB1(); // levels in ry unit
			levels_eV[i] = nanoparticles[i].getCB1()*Constants.rytoev; // levels in eV unit
		}
		Arrays.sort(levels);
		double[] levelCount = util.Utility.calcHistogram(levels, levels[0], levels[levels.length-1], bins);
		double[] levelLinear = util.Utility.linearSpread(levels[0], levels[levels.length-1], bins);
		
		//Fit a gaussian curve to the histogram
		GaussianCurveFitter fitter = GaussianCurveFitter.create();
	    WeightedObservedPoints obs = new WeightedObservedPoints();
	    for (int index = 0; index < levelLinear.length; index++) {
	            obs.add(levelLinear[index], levelCount[index]);}
	    double[] bestFit = fitter.fit(obs.toList());
	    // bestFit = norm, mean, sigma  ---> return sigma ---> FWHM
	    System.out.println("Sigma is: " + bestFit[2]*Constants.rytoev);
		
	    //Return the full width half max
	    FWHM = bestFit[2]*2*Math.sqrt(2*Math.log(2));
	    return FWHM; 
	}
	
	
	public void setupElectrons(){
		
		return;
	}
	

	public static void main(String[] args) {

		// Necessary parameters
		//(INT_t steps, INT_t sample_no, INT_t e_number, INT_t h_number, INT_t nanops_number, 
		//FLOAT_t large_ratio, FLOAT_t temp, FLOAT_t thr, str features, FLOAT_t voltage_ratio, output):
        Map<String, Object> params = new HashMap<>();
        
        params.put("nelec", 200);
    	params.put("nholes", 0);
    	params.put("nnanops", 400);
    	params.put("sample_no", 0);
    	params.put("feature", "mobility");
    	params.put("temperature", 80.0);
    	params.put("thr", 0.0);
    	params.put("np_diam", 5.0);
        
        Sample newsample = new Sample(params);
        
        
        
        
        
        //System.out.println(newsample.getElectronMass());
        //System.out.println(newsample.getFWHM()*Constants.rytoev);
        //Nanoparticle[] nanoparticles;
        //nanoparticles = setupNanoparticles(newsample);
        
        //System.out.println(nanoparticles[0].x/Constants.nmtobohr);
        //System.out.println(nanoparticles[399].getCB1());
		
		
		
	}

}
