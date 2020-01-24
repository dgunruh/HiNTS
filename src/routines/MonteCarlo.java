package routines;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import classes.Sample;
import util.Constants;
import util.MersenneTwisterFast;

public class MonteCarlo {

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
        MersenneTwisterFast r;

        r = new MersenneTwisterFast(new int[]{0x123, 0x234, 0x345, 0x456});
        
        for(int i=0; i<10; i++){
	        System.out.println(r.nextDouble());
	        
	        System.out.println(util.Constants.ry_ps);
	        
	        System.out.println(util.Configuration.hoppingMechanism);
	        
	        
	        Map<String, Object> params = new HashMap<>();
	        
	        params.put("nelec", 50);
	    	params.put("nholes", 0);
	    	params.put("nnanops", 400);
	    	params.put("sample_no", 0);
	    	params.put("feature", "mobility");
	    	params.put("temperature", 300.0);
	    	params.put("thr", 0.2);
	        
	        Sample newsample = new Sample(params);
	        
	        
	        System.out.println(newsample.getElectronMass());
        
        }
	}

}
