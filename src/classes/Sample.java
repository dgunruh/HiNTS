package classes;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import com.google.common.base.Charsets;
import com.google.common.io.Files;
import routines.Setup;
import util.Configuration;
import util.Constants;
import util.MersenneTwisterFast;

public class Sample {
    String capacitanceModel;
	//boolean lcapzunger, lcapdelerue; 
    String effectiveMedium;
    //boolean lemmg, lemla, lemlll, lempnone;
    String hoppingMechanism;
    //boolean lma, lmarcus; 
    String poissonSolver; 
    //boolean lpoissonnone, lpoissonewald, lpoissonnn;
    public String latticeStructure;
    public double latticeAngleGamma;
    public double latticeAngleAlpha;
    public double latticeAngleBeta;

    boolean perx, pery, perz, bimodal, twolayer, lcapacitance0, pennmodel = false;
    
    double marcusprefac2, jumpfreq, emass, hmass;
    double reorgenergy, bradii, capacitance,ligandlength, sourcewf, drainwf;
    double chemicalPotential = 0, V_prime;  // used for grand canonical ensemble

    public double temperature, voltage ;  
    public double near_n_dist_thr, nextn_n_dist_thr, neck_dist_thr, ediff_thr, ediff_neck_thr, beta, packingfraction;
    double ediff_hopping_thr; //for the overlap energy for different disorder values
    public int e_degeneracy, h_degeneracy, nbnd;
    public double cellx, celly, cellz, cellz_nm;
    public double rateOnSample = 0;
    public double systemEnergy = 0;
    long sampleCurrent = 0;
    public double elapsedTime = 0.0;
    
    public String folderName; //the folder where the nanoparticle data is stored
    public int nelec   ;   // number of electrons
    public int nelec_0 ;   // starting number of electrons
    public int nhole   ;   // number of holes
    public int expected_nnanops ;   // expected number of nanoparticles (eg 10x2x20 gives 400 nanops)
    public int nnanops; //actual number (if vacancies, grains, etc)
    int sample_number  ;   // store the nanoparticle file to open
    double closeNeighbor_thr; // percentage read from program
    double neckedNeighbor_thr; // percentage read from program
    public double screeningFactor; //how much to reduce charging energy by, if necessary
    public boolean randomSeed; //do we use a random seed for MersenneTwisterFast, or preset seed
    public boolean necking; //are the nanoparticles necked together
    public boolean reverse_the_sample; //should the sample characteristic be flipped around. Unnecessary if the sample was "pre-flipped" or if the sample is unflipped regardless
  
    double npdc, liganddc, metalPrefactor, FWHM;
    double average_spacing = 2*Configuration.ligandLength*Constants.nmtobohr;
    String filename, feature;
    
    //These are the parameters for the extra electron tracking
    public double totalZeroElectronNanoparticles;
    public double totalOneElectronNanoparticles;
    public double totalTwoElectronNanoparticles;
    public double totalThreeElectronNanoparticles;
    
    
    //
    Nanoparticle[] nanoparticles;
    ArrayList<Nanoparticle> sources;
    ArrayList<Nanoparticle> drains;
    HashSet<Nanoparticle> currentUpdate = new HashSet<Nanoparticle>();
    HoppingEvent currentEvent, previousEvent, latestEvent;
    public HashMap<Integer, Integer> nanoparticle_idList = new HashMap<Integer, Integer>();
    ArrayList<ArrayList<Object>> neighborList = new ArrayList<ArrayList<Object>>();
    
    public double leftmost_NP_z = 36.8*Constants.nmtobohr, rightmost_NP_z = 126.4*Constants.nmtobohr;
    
    MersenneTwisterFast rng;

 
    public Sample(Map<String,Object> params){
    	
    	nelec = (int) params.get("nelec");
    	nelec_0 = nelec;
    	nhole = (int) params.get("nholes");
    	expected_nnanops = (int) params.get("expected_nnanops");
    	nnanops = expected_nnanops; //default value
    	sample_number = (int) params.get("sample_no");
    	reverse_the_sample = (boolean) params.get("reverse_the_sample");
    	feature = (String) params.get("feature");
    	temperature = (double) params.get("temperature") * Constants.kelvintory;
    	closeNeighbor_thr = (double) params.get("closeNeighbor_thr");
    	neckedNeighbor_thr = (double) params.get("neckedNeighbor_thr");
    	screeningFactor = (double) params.get("screeningFactor");
    	randomSeed = (boolean) params.get("randomSeed");
    	necking = (boolean) params.get("necking");
    	folderName = (String) params.get("folderName");
    	
    	if(feature=="mobility"){
    		if(!reverse_the_sample) {
    			//voltage = 5*30*Constants.kelvintory*0.1*20 / Constants.sqrt2; //Multiplied by 8. This was for the moule samples
    			voltage = 5*30*Constants.kelvintory*0.1*20 / Constants.sqrt2;
    			System.out.println("Mobility run, voltage is "+voltage*Constants.rytoev);
    		}
    		else {
    			voltage = -5*30*Constants.kelvintory*0.1*20 / Constants.sqrt2;
    			System.out.println("Negative voltage");
    		}
    		//voltage = 0.5*30*Constants.kelvintory*0.1*25 / Constants.sqrt2; //25 nanoparticles, so each get this 30K voltage
    	}
    	if(feature=="iv")
    		voltage = 30*Constants.kelvintory*((double)params.get("voltage_ratio"))*.1*20/Constants.sqrt2;
    	
    	

    	loadConfiguration();        
   
		nanoparticles = Setup.setupNanoparticles(this);
		
		/*
		if(Configuration.latticeStructure == "Moule") {
			redefineSourcesAndDrains();
		}
		*/
		
		System.out.println("Number of Sources: " + sources.size());
		System.out.println("Number of Drains: " + drains.size());
		if(folderName.substring(folderName.length() - 4).equals(".00D")) {
			FWHM = 0.0;
		}
		else {
			FWHM = Setup.getFWHM(nanoparticles, this);
		}
        System.out.println("FWHM is: " + FWHM*Constants.rytoev);
		//FWHM = 0.1;
        ediff_thr = FWHM * closeNeighbor_thr;
        //ediff_neck_thr = FWHM * neckedNeighbor_thr;
        ediff_neck_thr = .014 * Constants.evtory;
        //ediff_hopping_thr = .014*Constants.evtory; //an overlap energy of 14 meV
        ediff_hopping_thr = .014*Constants.evtory;
        if(necking) get_necked_neighbors();
        //get_necked_neighbors();
        buildNeighborList(false, true);
        dielectrics(true);
        set_selfenergy();
        
        // throw electrons on to nanoparticles
        for(int i=0; i<nelec; i++)
        	throw_electron();
        
        //This method fills every NP, then adds on to the even landscape randomly
//        for(Nanoparticle np: nanoparticles) {
//        	Electron e = new Electron();
//        	e.setHost(np, 0);
//        	np.add_electron(e, 0);
//        	double[] destination = {np.x, np.z};
//        	e.visitedNPs.add(destination);
//        }
//        for(int i = 0; i < nelec%800; i++) throw_electron();
 
        //throw holes on to nanoparticles
        for(int i=0; i < nhole; i++)
        	throw_hole();
        
        //simulation();
        
        
        
    }
    
    private void redefineSourcesAndDrains() {
    	sources.clear();
    	drains.clear();
    	double tempMax = 0;
    	double tempMin = 10000;
    	for(Nanoparticle NP: nanoparticles) {
    		double NP_z = NP.z;
    		if (NP_z < tempMin) {
    			tempMin = NP_z; 
    		}
    		else if (NP_z > tempMax) {
    			tempMax = NP_z;
    		}
    		NP.source = false;
    		NP.drain = false;
    	}
    	
    	for(Nanoparticle NP: nanoparticles) {
    		double z = NP.z, radius = NP.radius;
    		if(z-radius <= (tempMin + near_n_dist_thr)){
				NP.source = true;
				sources.add(NP);
			}
			if(z+radius >= (tempMax - near_n_dist_thr)){
				NP.drain = true;
				drains.add(NP);
			}
    	}
    	
    }
    
    private void loadConfiguration()  {
    	
    	// setting up the random number generator
    	if(randomSeed) {
    		int length = 20;
    		int rng_seed = (int)(Math.random()*length);
    		rng = new MersenneTwisterFast(rng_seed);
    	}
    	else {
    		rng = new MersenneTwisterFast(Configuration.rngSeed);
    	}
    	
    	emass = Configuration.emass;
    	hmass = Configuration.hmass;
    	ligandlength = Configuration.ligandLength*Constants.nmtobohr;
    	capacitance = Configuration.capacitance*2.0/Constants.evtory; // used for zunger
    	near_n_dist_thr = Configuration.near_n_distThr*Constants.nmtobohr;
    	nextn_n_dist_thr = Configuration.nextn_n_distThr*Constants.nmtobohr;
    	neck_dist_thr = Configuration.neckDistThr*Constants.nmtobohr;
    	liganddc = Configuration.ligandDC;
    	npdc = Configuration.npDC + (temperature*Constants.ry_to_kelvin - 300.0)*Configuration.epsilon_per_K;
    	nbnd = Configuration.nBands;
    	e_degeneracy = Configuration.e_degeneracy;
    	h_degeneracy = Configuration.h_degeneracy;
    	
    	metalPrefactor = Configuration.metallicPrefactor;
    	
    	perx = Configuration.PERX;
    	pery = Configuration.PERY;
    	perz = Configuration.PERZ;
    	
    	latticeStructure = Configuration.latticeStructure;
    	latticeAngleGamma = Configuration.latticeAngleGamma;
    	latticeAngleAlpha = Configuration.latticeAngleAlpha;
    	latticeAngleBeta = Configuration.latticeAngleBeta;
    	sources = new ArrayList<Nanoparticle>();
    	drains = new ArrayList<Nanoparticle>();
    	
    	hoppingMechanism = Configuration.hoppingMechanism;
    	effectiveMedium = Configuration.effectiveMedium;
    	capacitanceModel = Configuration.capacitanceModel;
    	poissonSolver = Configuration.poissonSolver;
    	
    	lcapacitance0 = Configuration.lcapacitance0;
    	marcusprefac2 = Configuration.marcusprefac2*Constants.evtory*Constants.evtory;
    	jumpfreq = Configuration.jumpFreq*Constants.ry_ps;
    	
    	reorgenergy = Configuration.reorgenergy*Constants.evtory;
    	bradii = Configuration.bradii*Constants.nmtobohr;
    	
    	bimodal = Configuration.biModal;
    	twolayer = Configuration.twoLayer;
    	
    	if(!bimodal){
    		//String prefix = "c:/Users/Davis/Research/Superlattice Project/npSolids/";
    		String prefix = "c:/Users/dunru/GitHub/HiNTS/data/";
    		String middle = folderName; //_.03D"; // looking at grain boundaries
    		String end = "/npSample" + sample_number + ".inp";
    		if(necking) {
    			end = "/npSample" + sample_number + ".inp";
    		}
    		//String end = "/npSample" + sample_number + ".inp";
    		
    		/* Moule project
    		String prefix= "./data/14x14x6/";
    		//String middle = String.valueOf(nnanops)+"_"+Configuration.diameter;
    		String middle = "npSample";
    		String end = String.valueOf(sample_number)+ ".inp";
    		filename = prefix + middle + end;
    		
    		if(Configuration.latticeStructure == "Moule") {
    			
	    		prefix= "./data/";
	    		//String middle = String.valueOf(nnanops)+"_"+Configuration.diameter;
	    		middle = "QD_com_xyz_1846_trimmed_square3";
	    		if(sample_number == 1) {
	    			middle = "QD_com_xyz_1846_trimmed_square3_inverted";
	    		}
	    		end = ".txt";
    		}
    		*/
    		
    		filename = prefix + middle + end;
    		System.out.println(filename);
    	}
    	
    	//get actual number of nanoparticles
    	try
    	(
    	   FileReader       input = new FileReader(filename);
    	   LineNumberReader count = new LineNumberReader(input);
    	)
    	{
    	   while (count.skip(Long.MAX_VALUE) > 0)
    	   {
    	      // Loop just in case the file is > Long.MAX_VALUE or skip() decides to not read the entire file
    	   }
    	   
    	   nnanops = count.getLineNumber() - 7; // +1 because line index starts at 0, -7 for starting lines, -1 because ends in empty line
    	   if(necking) {
    		   nnanops -= 1;
    	   }
    	} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
    	System.out.println("Number of nanoparticles is: " + nnanops);
    	
		List<String> lines;
		try {
			lines = Files.readLines(new File(filename), Charsets.UTF_8);
			
			/* For Moule style samples
			cellx = Double.valueOf(Arrays.asList(lines.get(0).split(",")).get(1))*Constants.nmtobohr;
			celly = Double.valueOf(Arrays.asList(lines.get(1).split(",")).get(1))*Constants.nmtobohr;
			cellz = Double.valueOf(Arrays.asList(lines.get(2).split(",")).get(1))*Constants.nmtobohr;
			if(Configuration.latticeStructure == "Moule") {
				celly = Double.valueOf(Arrays.asList(lines.get(0).split(",")).get(1))*Constants.nmtobohr;
				cellz = Double.valueOf(Arrays.asList(lines.get(1).split(",")).get(1))*Constants.nmtobohr;
				cellx = Double.valueOf(Arrays.asList(lines.get(2).split(",")).get(1))*Constants.nmtobohr;
			}
			*/
			cellx = Double.valueOf(Arrays.asList(lines.get(3).split(",")).get(1))*Constants.nmtobohr;
			celly = Double.valueOf(Arrays.asList(lines.get(4).split(",")).get(1))*Constants.nmtobohr;
			cellz = Double.valueOf(Arrays.asList(lines.get(5).split(",")).get(1))*Constants.nmtobohr;
			cellz_nm = cellz/Constants.nmtobohr;
			
			//leftmost_NP_z = cellz;
			
			previousEvent = new HoppingEvent();
			currentEvent = new HoppingEvent();
			latestEvent = new HoppingEvent();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// calculate V_prime and chemical potential for Grand Canonical
		//V_prime = (cellx*celly*cellz) / (Math.pow(2*Math.PI, 1.5) / Math.pow(Configuration.emass*Constants.k_boltzmann_ry*temperature, -1.5)); 
		double lambda = Math.pow(2*Math.PI/(Configuration.emass*temperature), 0.5); //lambda in bohr (bohrs are 1 in Ry units)
		System.out.println("Lambda is: " + lambda);
		V_prime = (cellx*celly*cellz) / Math.pow(lambda, 3);
		chemicalPotential = temperature*Math.log(nelec/V_prime); //chemical potential of ideal gas of electrons
		//chemicalPotential = -0.5*Constants.evtory;
		System.out.println("V_prime is "+ V_prime+" "+cellx);
		//System.out.println("Chemical potential is: " + chemicalPotential);

	}
    
    public double getElectronMass(){
    	return emass;
    }
    
    public String getFilename(){
    	return filename;
    }
    
    public double getFWHM() {
    	return FWHM;
		
	}

    private void dielectrics(boolean connected_z) {
		
    	//Determine packing fraction and average NP diameter
    	double npVolume = 0;
    	double average_diameter = 0;
    	for(int i=0; i<nanoparticles.length; i++){
    		npVolume += 4.0/3.0*Constants.pi*Math.pow(nanoparticles[i].radius, 3);
    		average_diameter += nanoparticles[i].diameter;
    		if(pennmodel == true)
    			nanoparticles[i].dcin = 1+(npdc-1)/(1+Math.pow((bradii/nanoparticles[i].diameter),2));
    		else
    			nanoparticles[i].dcin = npdc;
    	}
    	average_diameter = average_diameter/nanoparticles.length;
    	packingfraction = npVolume / (cellx*celly*cellz);
    	
    	// set dcout for each NP individually
    	switch (effectiveMedium){
			case "mg":
				for(int i=0; i<nanoparticles.length; i++) {
					nanoparticles[i].dcout = liganddc*(nanoparticles[i].dcin*(1+2*packingfraction)-liganddc*(2*packingfraction-2))/(liganddc*(2+packingfraction)+nanoparticles[i].dcin*(1-packingfraction));
				}
				
				//dcout = liganddc*(npdc*(1+2*packingfraction)-liganddc*(2*packingfraction-2))/(liganddc*(2+packingfraction)+npdc*(1-packingfraction));
				//dcout = sample.liganddc*(sample.npdc*(1+2*packingfraction)-sample.liganddc*(2*packingfraction-2))/(sample.liganddc*(2+packingfraction)+sample.npdc*(1-packingfraction))
				break;
			
			case "shklovskii global":
				double delta = average_diameter*Math.pow(Math.PI,6.0/5.0)/Math.pow(2, 11.0/5.0)*Math.pow(liganddc/npdc,6.0/5.0);
				if(average_spacing >=0) {
					for(int i=0; i<nanoparticles.length; i++) {
						nanoparticles[i].dcout = Math.PI/2*liganddc*Math.pow(average_diameter/(2*average_spacing + 2*delta), 1.0/3.0);
					}
				}
				else {
					for(int i=0; i<nanoparticles.length; i++) {
						nanoparticles[i].dcout = npdc*Math.sqrt(2*(Math.abs(average_spacing)+delta)/average_diameter);
					}
				}
				break;

			case "shklovskii local":
		    	//Average the inter-NP spacing over next-nearest-neighbors and nearest-neighbors
				//This works for non-necked NPs, or necked NPs where the overlap is equal to the neck diameter
		    	for(int i=0; i<nanoparticles.length; i++){
		        	double totalSpacing = 0.0;
		        	double totalDiameter = 0.0;
		        	int totalSpaces = 0;
		        	int visitedNPs = 0;
		    		ArrayList<Nanoparticle> visitedNeighbors = new ArrayList<Nanoparticle>();
		    		double totalVolume = 4.0/3.0*Constants.pi*Math.pow(nanoparticles[i].radius, 3);
		    		
		    		for(Nanoparticle neighborNP : nanoparticles[i].nextNearestNeighbors) {
		    			visitedNeighbors.add(neighborNP);
		    			totalDiameter += neighborNP.diameter;
		    			visitedNPs += 1;
		    			
		    			//calculate average spacing directly
		    			if(!necking) {
			    			for ( Map.Entry<Nanoparticle, Double> entry : neighborNP.edgeDistanceMap.entrySet()) {
								Nanoparticle neighbor = entry.getKey();
								Double spacing = entry.getValue();
								if (!visitedNeighbors.contains(neighbor)) {
									totalSpacing += spacing;
									totalSpaces += 1;
								}
							}
		    			}
		    			
		    			//calculate average spacing from neck widths (assume two overlapping spheres)
		    			boolean simplify = true;
		    			if(necking && simplify) {
		    				for(Nanoparticle newNeighborNP : neighborNP.nearestNeighbors) {
		    					if(!visitedNeighbors.contains(newNeighborNP)) {
		    						//Check if there is a neck
		    						//If neck, then calculate artificial spacing
		    						if(neighborNP.neckRadiusMap.containsKey(newNeighborNP)) {
		    							double w = neighborNP.neckRadiusMap.get(newNeighborNP);
		    							double R = neighborNP.radius;
		    							double r = newNeighborNP.radius;
		    							double spacing = Math.sqrt(Math.pow(R, 2.0) - Math.pow(w, 2.0)) 
		    									+ Math.sqrt(Math.pow(R, 2.0) + Math.pow(r, 2.0) - 2*Math.pow(w, 2.0));
		    							totalSpacing += spacing;
		    							totalSpaces += 1;
		    						}
		    						//If no neck, then use actual spacing
		    						else {
		    							double spacing = neighborNP.edgeDistanceMap.get(newNeighborNP);
		    							totalSpacing += spacing;
		    							totalSpaces += 1;
		    						}
		    					}
		    				}
		    			}
		    			
		    			//calculate filling fraction, which can then be translated to an average spacing number
		    			else if(necking && !simplify) {
		    				totalVolume += 4.0/3.0*Constants.pi*Math.pow(neighborNP.radius, 3);
		    				for(Nanoparticle newNeighborNP : neighborNP.nearestNeighbors) {
		    					if(!visitedNeighbors.contains(newNeighborNP)) {
		    						if(neighborNP.neckRadiusMap.containsKey(newNeighborNP)) {
		    							//get the neck width. Assume this is the correct neck width no matter what the actual
		    							//NP-NP center-to-center spacing is
		    							double w = neighborNP.neckRadiusMap.get(newNeighborNP);
		    							
		    							//Now, we need the length of the cylinder that we will represent the neck as
		    							double r_1 = neighborNP.radius;
		    							double r_2 = newNeighborNP.radius;
		    							double d_1 = Math.sqrt(Math.pow(r_1, 2.0) - Math.pow(w, 2.0));
		    							double d_2 = Math.sqrt(Math.pow(r_2, 2.0) - Math.pow(w, 2.0));
		    							double w_length = neighborNP.centerDistanceMap.get(newNeighborNP) - d_1 - d_2;
		    							
		    							//Now we get the volume of the intersecting spherical caps
		    							double cap_1 = r_1 - d_1;
		    							double cap_2 = r_2 - d_2;
		    							double cap_1_V = Constants.pi*Math.pow(cap_1, 2.0)*(3*r_1-cap_1);
		    							double cap_2_V = Constants.pi*Math.pow(cap_2, 2.0)*(3*r_2-cap_2);
		    							
		    							//Finally, get the additional volume that the neck provides
		    							totalVolume += (Constants.pi*Math.pow(w, 2.0)*w_length - cap_1_V - cap_2_V);
		    						}
		    					}
		    				}
		    			}
		    		}
					double averageSpacing = totalSpacing/totalSpaces;
					double averageDiameter = totalDiameter/visitedNPs;
					double del = averageDiameter*Math.pow(Math.PI,6.0/5.0)/Math.pow(2, 11.0/5.0)*Math.pow(liganddc/nanoparticles[i].dcin,6.0/5.0);
					if(averageSpacing >=0) {
						nanoparticles[i].dcout = Math.PI/2*liganddc*Math.pow(averageDiameter/(2*averageSpacing + 2*del), 1.0/3.0);
					}
					else {
						nanoparticles[i].dcout = nanoparticles[i].dcin*Math.sqrt(2*(Math.abs(averageSpacing)+del)/averageDiameter);
					}

		    	}
		    	break;
		    	
			case "shklovskii sphere":

		    	//Compute the spatial filling of the local area around each NP.
		    	//To do this, create a spherical bounding box, and calculate the volume overlap
		    	//of the NPs and necks with this box. Note: necks are approximated as spheres for simplicity. 
		    	double f, j, x, s, d;
		    	for(int i=0; i<nanoparticles.length; i++){
		    		f = calculate_local_filling_fraction(nanoparticles[i], 2*Configuration.nextn_n_distThr*Constants.nmtobohr);
		    		//System.out.println("Filling fraction: " + f);
		    		//Now turn the filling fraction into Shklovskii's s/d ratio
		    		x = -0.4;
		    		j = 1.0;
		    		while (x < 0.4 && !(Math.abs(j) < .0001)) {
		    			j = -6.2832 + 12*f + 36*f*x + (28.27 + 36*f)*Math.pow(x, 2.0) + (9.4 + 12*f)*Math.pow(x, 3.0);
		    			x += .00001;
		    		}
		    		//System.out.println("s/d ratio: " + x);
		    		//Two options to from s/d to s and d separately. 1: use NP d. 2: use average NP d in bounding box. Here I use (1)
		    		s = x*nanoparticles[i].diameter;
		    		d = nanoparticles[i].diameter;
		    		delta = d*Math.pow(Math.PI,6.0/5.0)/Math.pow(2, 11.0/5.0)*Math.pow(liganddc/npdc,6.0/5.0);
					if(s >=0) {
						nanoparticles[i].dcout = Math.PI/2*liganddc*Math.pow(d/(2*s + 2*delta), 1.0/3.0);
					}
					else {
						nanoparticles[i].dcout = npdc*Math.sqrt(2*(Math.abs(s)+delta)/d);
					}
					//System.out.println("Dielectric constant: " + nanoparticles[i].dcout);
		    	}
		    	
		    	break;
				
			default:
				for(int i=0; i<nanoparticles.length; i++) {
					nanoparticles[i].dcout = liganddc;
				}
				System.out.println("Entered default");
				break;
		}	
	}
    
    //Calculate volumetric filling fraction of a NP in a spherical bounding box.
    //Inputs: bounding radius of the spherical box
    private double calculate_local_filling_fraction(Nanoparticle np, double boundingRadius) {
    	//We want to add up the entire volume of nanoparticles (and necks if applicable) inside a spherical box
    	//Using a cuboid as the box is much more difficult, and not necessary!
    	double totalVolume = 4.0*Constants.pi*Math.pow(np.radius, 3.0)/3.0;
    	ArrayList<Nanoparticle> visitedNanoparticles = new ArrayList<Nanoparticle>();
    	for(int i=0; i<nanoparticles.length; i++){
    		if(nanoparticles[i]!=np) {
    			visitedNanoparticles.add(nanoparticles[i]);
    			double r = nanoparticles[i].radius;
    			double R = boundingRadius;
    			double d = np.npnpdistance(nanoparticles[i], this, true, true);
    			if(d == 0) {System.out.println("zero d!");}
    			double overlapVolume = volumeSphericalOverlap(r, R, d);
    			totalVolume += overlapVolume;
    			boolean inside = overlapVolume > 0.0;
    			//check if it connects to any nanoparticles inside the bounding sphere
    			if(necking) {
    				for ( Map.Entry<Nanoparticle, Double> entry : nanoparticles[i].neckRadiusMap.entrySet()) {
						Nanoparticle neighbor = entry.getKey();
						if(!visitedNanoparticles.contains(neighbor)) {
							double neckRadius = entry.getValue();
							double dneighbor = np.npnpdistance(neighbor, this, true, true);
							double equivalentRadius = 0.0; //the radius of the "neck sphere"
							
							//if NP_i lives inside the bounding sphere, will check neck regardless
							if (inside || dneighbor - neighbor.radius <= R) {
								//Get volume of non-overlapping portion of neck
								double d_1 = Math.sqrt(Math.pow(r, 2.0) - Math.pow(neckRadius, 2.0));
								double d_2 = Math.sqrt(Math.pow(neighbor.radius, 2.0) - Math.pow(neckRadius, 2.0));
								double w_length = nanoparticles[i].centerDistanceMap.get(neighbor) - d_1 - d_2;
								
								//Now we get the volume of the intersecting spherical caps
								double cap_1 = r - d_1;
								double cap_2 = neighbor.radius - d_2;
								double cap_1_V = Constants.pi*Math.pow(cap_1, 2.0)*(3*r-cap_1);
								double cap_2_V = Constants.pi*Math.pow(cap_2, 2.0)*(3*neighbor.radius-cap_2);
								
								//Finally, get the additional volume that the neck provides
								double neckVolume = (Constants.pi*Math.pow(neckRadius, 2.0)*w_length - cap_1_V - cap_2_V);
								equivalentRadius = Math.pow(3.0*neckVolume/(4*Constants.pi), 1.0/3.0);
								
							}
							
							if (equivalentRadius > 0.0) {
								//get coordinates for neck center
								double deltax = neighbor.x - nanoparticles[i].x;
								if (Math.abs(deltax) > Math.abs(neighbor.x-(nanoparticles[i].x-cellx)))
									deltax = neighbor.x-(nanoparticles[i].x-cellx);
								if (Math.abs(deltax) > Math.abs(neighbor.x-(nanoparticles[i].x+cellx)))
									deltax = neighbor.x-(nanoparticles[i].x+cellx);
								
								double deltay = neighbor.y - nanoparticles[i].y;
								if (Math.abs(deltay) > Math.abs(neighbor.y-(nanoparticles[i].y-celly)))
									deltay = neighbor.y-(nanoparticles[i].y-celly);
								if (Math.abs(deltay) > Math.abs(neighbor.y-(nanoparticles[i].y+celly)))
									deltay = neighbor.y-(nanoparticles[i].y+celly);
								
								double deltaz = neighbor.z - nanoparticles[i].z;
								if (Math.abs(deltaz) > Math.abs(neighbor.z-(nanoparticles[i].z-cellz)))
									deltaz = neighbor.z-(nanoparticles[i].z-cellz);
								if (Math.abs(deltaz) > Math.abs(neighbor.z-(nanoparticles[i].z+cellz)))
									deltaz = neighbor.z-(nanoparticles[i].z+cellz);
								double[] neckCenter = new double[]{nanoparticles[i].x + deltax/2.0, nanoparticles[i].y + deltay/2.0, nanoparticles[i].z + deltaz/2.0};
								double neckDistance = np.nppointdistance(neckCenter, this);
								double neckOverlapVolume = volumeSphericalOverlap(equivalentRadius, r, neckDistance);
								totalVolume += neckOverlapVolume;
							}
						}
    				}
    			}
    		}
    		
    	}
    	
    	return totalVolume*3.0/(4.0*Constants.pi*Math.pow(boundingRadius, 3.0));
    }
    
    //Inputs: 
    //r: the radius of the smaller sphere overlapping the bounding box
    //R: the radius of the larger sphere (the bounding box)
    //d: the distance between the centers of the two spheres
    public double volumeSphericalOverlap(double r, double R, double d) {
    	double volume = 0.0;
		//check if the nanoparticle lives fully inside bounding sphere
		if(d + r <= R) {
			volume = 4.0*Constants.pi*Math.pow(r, 3.0)/3.0;
		}
		//check if the nanoparticles lives partially inside the bounding sphere, and center outside sphere
		else if(d - r <= R && d > R) {
			double h1 = R - (Math.pow(R, 2.0) - Math.pow(r, 2.0) + Math.pow(d, 2.0))/(2*d);
			double h2 = (Math.pow(R, 2.0) - Math.pow(r, 2.0) + Math.pow(d, 2.0))/(2*d) - (d-r);
			volume = Constants.pi*Math.pow(h1, 2.0)*(3*R - h1)/3.0 + Constants.pi*Math.pow(h2, 2.0)*(3*r - h2)/3.0;
		}
		//check if the nanoparticle lives partially inside the bounding sphere and center inside sphere
		else if(d - r <= R && d < R) {
			double h1 = R - (Math.pow(R, 2.0) - Math.pow(r, 2.0) + Math.pow(d, 2.0))/(2*d);
			double h2 = r - (Math.pow(R, 2.0) - Math.pow(r, 2.0) - Math.pow(d, 2.0))/(2*d);
			double volumeOutside = Constants.pi*Math.pow(h2, 2.0)*(3*r - h2)/3.0 - Constants.pi*Math.pow(h1, 2.0)*(3*R - h1)/3.0;
			volume = 4.0*Constants.pi*Math.pow(r, 3.0)/3.0 - volumeOutside;
		}
		
		return volume;
    }

    //works for one band only
    private void set_selfenergy(){
		for(int i=0; i<nanoparticles.length; i++){
			
			switch (capacitanceModel) {
				case "delerue":
					//System.out.println(i);
					nanoparticles[i].selfenergy0 = ((Constants.e2/nanoparticles[i].radius)*((0.5/nanoparticles[i].dcout - 0.5/nanoparticles[i].dcin) + 0.47*(nanoparticles[i].dcin-nanoparticles[i].dcout)/((nanoparticles[i].dcin+nanoparticles[i].dcout)*nanoparticles[i].dcin))); 
					//System.out.println("sigma0 = " + nanoparticles[i].selfenergy0);
					nanoparticles[i].selfenergy = (Constants.e2/nanoparticles[i].radius*(1.0 / nanoparticles[i].dcout + 0.79 /nanoparticles[i].dcin));
					
					//System.out.println("sigma = " + nanoparticles[i].selfenergy);
					
					//make hole energies negative
					//nanoparticles[i].hselfenergy0 = -(Constants.e2/nanoparticles[i].radius*((0.5 /nanoparticles[i].dcout-0.5 /nanoparticles[i].dcin)+0.47 /nanoparticles[i].dcin*(nanoparticles[i].dcin - nanoparticles[i].dcout)/(nanoparticles[i].dcin + nanoparticles[i].dcout)));
					//nanoparticles[i].hselfenergy = -(Constants.e2/nanoparticles[i].radius*(1.0 / nanoparticles[i].dcout + 0.79 /nanoparticles[i].dcin));
					nanoparticles[i].hselfenergy0 = -nanoparticles[i].selfenergy0;
					nanoparticles[i].hselfenergy = -nanoparticles[i].selfenergy;
				break;
				

			case "zunger":
					nanoparticles[i].selfenergy0 = Constants.e2/(2.0 *capacitance*nanoparticles[i].radius);
					nanoparticles[i].selfenergy = Constants.e2/(2.0 * capacitance*nanoparticles[i].radius);
					
					//make hole energies negative
					nanoparticles[i].hselfenergy0 = -Constants.e2/(2.0 *capacitance*nanoparticles[i].radius);
					nanoparticles[i].hselfenergy= -Constants.e2/(2.0 *capacitance*nanoparticles[i].radius); //charging energy
				break;
			}
		}
	}
    
    public void get_necked_neighbors() {
    	String prefix= "c:/Users/Davis/Research/Superlattice Project/neckDistributions/" + folderName;
		String middle = "/neckSample";
		String end = String.valueOf(sample_number/2)+ ".inp";
		
		if (Configuration.latticeStructure == "Moule") {
			prefix = "./data/";
	    	middle = "QD_necking_1846_trimmed_square";
	    	end = ".txt";
		}
		//System.out.println("Neck sample number: " + (sample_number+1000)/2);
		String neck_filename = prefix + middle + end;
		
		for(int i = 0; i < nanoparticle_idList.size(); i++) {
			neighborList.add(new ArrayList<Object>());
		}
    	
    	//ArrayList<Object> nn_neck = new ArrayList<Object>();
    	int current_id = 1;
    	
    	
    	//neighborList will contain arrays which list the ids (and corresponding necks) of necked neighbors. These ids will be
    	//the order positions in the nanoparticles[] array
    	//the order of neighborList will match nanoparticles[]
    	List<String> lines;
    	try {
			lines = Files.readLines(new File(neck_filename), Charsets.UTF_8);
			int i = 0;
			for(String lineWhole : lines) {
				if(i > 0) {
					List<String> line = Arrays.asList(lineWhole.split(" "));
					int thisID = Integer.valueOf(line.get(0));
					int neighborID = Integer.valueOf(line.get(1));
					double neck_diameter = Double.valueOf(line.get(3));
					
					int neighborList_neighborID = 1;
					try {
						neighborList_neighborID = nanoparticle_idList.get(neighborID);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						System.out.println("Neighbor id: " + neighborID);
						System.out.println(nanoparticle_idList);
					}
					int neighborList_thisID = nanoparticle_idList.get(thisID); 
					ArrayList<Object> nn_necks = neighborList.get(neighborList_thisID);
					ArrayList<Object> neighbor_nn_necks = neighborList.get(neighborList_neighborID);
					
					nn_necks.add(neighborList_neighborID);
					nn_necks.add(neck_diameter);
					
					neighbor_nn_necks.add(neighborList_thisID);
					neighbor_nn_necks.add(neck_diameter);
					
					if(thisID != current_id) {
						//int current_neighborList_id = nanoparticle_idList.get(current_id); 
						//neighborList.set(current_neighborList_id, nn_neck);
						//current_id = this_id;
						
					}
				}
				i++;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    /**
	@param direction_selection, periodic in z
	@return nothing, modifies the sample neighborlist
	@throws what kind of exception does this method throw
	*/
    private void buildNeighborList(boolean directionSelection, boolean connected_z) {
    	int totalnn, totalcn, totalhcn, total_neckedNeighbors;
    	double edgeDist, centerDist;
    	ArrayList<Nanoparticle> nearestNeighbors, nextNearestNeighbors, closeNeighbors, closeHoleNeighbors, neckedNeighbors;
    	
    	totalcn=0;
    	totalnn=0;
    	totalhcn=0;
    	total_neckedNeighbors = 0;
    	double total_spacing = 0;
    	int total_spaces = 0;
    	
		for(int i=0; i<nanoparticles.length; i++){
			//total_neckedNeighbors = 0;
			
			nearestNeighbors = new ArrayList<Nanoparticle>();
			nextNearestNeighbors = new ArrayList<Nanoparticle>();
			closeNeighbors = new ArrayList<Nanoparticle>();
			closeHoleNeighbors = new ArrayList<Nanoparticle>();
			neckedNeighbors = new ArrayList<Nanoparticle>();
			nanoparticles[i].edgeDistanceMap = new HashMap<Nanoparticle, Double>();
			nanoparticles[i].centerDistanceMap = new HashMap<Nanoparticle, Double>();
			nanoparticles[i].electronHopMap = new HashMap<Nanoparticle, Integer>();
			nanoparticles[i].neckRadiusMap = new HashMap<Nanoparticle, Double>();
			
			for(int j=0; j<nanoparticles.length; j++){
				if(i != j){
					edgeDist = nanoparticles[i].npnpdistance(nanoparticles[j], this, connected_z, false);
					//if(edgeDist <= 0) {System.out.println("zero or negative edgeDist!!: " + edgeDist*Constants.bohrtonm);}
					centerDist = nanoparticles[i].npnpdistance(nanoparticles[j], this, connected_z, true);
					//if(centerDist <= 0) {System.out.println("zero or negative centerDist!!");}
					
					// Add nearest neighbors
					if(edgeDist <= near_n_dist_thr){
						total_spacing += edgeDist;
						total_spaces += 1;
						//if(edgeDist >=0) {
						nanoparticles[i].edgeDistanceMap.put(nanoparticles[j], edgeDist);
						nanoparticles[i].centerDistanceMap.put(nanoparticles[j], centerDist);
						//nanoparticles[i].electronHopMap.put(nanoparticles[j], 0);
						nearestNeighbors.add(nanoparticles[j]);
						//totalnn += 1;
						//add electron close neighbors
						if(Math.abs(nanoparticles[i].cbenergy[0]-nanoparticles[j].cbenergy[0]) < ediff_thr){
							closeNeighbors.add(nanoparticles[j]);
							totalcn += 1;
						}
						//add hole close neighbors
						if(Math.abs(nanoparticles[i].vbenergy[0]-nanoparticles[j].vbenergy[0]) < ediff_thr){
							closeHoleNeighbors.add(nanoparticles[j]);
							totalhcn += 1;
						}
						//}
						
					}
					
					if(edgeDist <= nextn_n_dist_thr) {
						nextNearestNeighbors.add(nanoparticles[j]);
					}
				}
			}
			
			if(necking) {
				ArrayList<Object> thisNeckNeighborList = neighborList.get(i);
				if(!thisNeckNeighborList.isEmpty()) {
					for(int j=0; j<thisNeckNeighborList.size()/2; j++){
						int iterator = 2*j;
						int index = (Integer) thisNeckNeighborList.get(iterator);
						//System.out.println("Nanoparticle id is: " + nanoparticles[i].id + " and it is necked to: " + nanoparticles[index].id);
						//int k = nanoparticle_idList.get(index);
						double neckDiameter = (Double) thisNeckNeighborList.get(iterator+1);
						double neckRadius = neckDiameter*Constants.nmtobohr/2.0;
						
						nanoparticles[i].neckRadiusMap.put(nanoparticles[index], neckRadius);
						//nanoparticles[k].neckRadiusMap.put(nanoparticles[i], neckRadius);
						neckedNeighbors.add(nanoparticles[index]);
						nearestNeighbors.remove(nanoparticles[index]); //no hopping transport when necks exist
						if(!nanoparticles[i].edgeDistanceMap.containsKey(nanoparticles[index])) {
							nanoparticles[i].edgeDistanceMap.put(nanoparticles[index], nanoparticles[i].npnpdistance(nanoparticles[index], this, connected_z, false));
						}
						if(!nanoparticles[i].centerDistanceMap.containsKey(nanoparticles[index])) {
							nanoparticles[i].centerDistanceMap.put(nanoparticles[index], nanoparticles[i].npnpdistance(nanoparticles[index], this, connected_z, true));
						}
						//nanoparticles[i].edgeDistanceMap.remove(nanoparticles[index]);
						//nanoparticles[i].centerDistanceMap.remove(nanoparticles[index]);
						//total_neckedNeighbors += 1;
					}
					
					//System.out.println("Nanoparticle id is: " + nanoparticles[i].id + " and the number of necked neighbors is: " + thisNeckNeighborList.size()/2);
				}
			}
			/*
			if(!nearestNeighbors.isEmpty()) {
				for(Nanoparticle neighborNP : new ArrayList<Nanoparticle>(nearestNeighbors)) {
					double edgeEdgeDist = nanoparticles[i].edgeDistanceMap.get(neighborNP);
					if(edgeEdgeDist < 0) {
						neckedNeighbors.add(neighborNP);
						int neighborIndex = nanoparticle_idList.get(neighborNP.id);
						double d = nanoparticles[i].centerDistanceMap.get(neighborNP);
						double R = nanoparticles[i].diameter/2;
						double r = nanoparticles[neighborIndex].diameter/2;
						//System.out.println("d is: " + d + "R is: " +R + " r is: " + r);
						double neckWidth = Math.sqrt(4*Math.pow(d*R, 2.0)-Math.pow(Math.pow(d, 2.0) + Math.pow(R, 2.0) - Math.pow(r, 2.0), 2.0))/(2*d);
						//System.out.println("neckWidth is: " + neckWidth);
						//double neckWidth = 0;
						nanoparticles[i].neckRadiusMap.put(neighborNP, neckWidth);
						nanoparticles[i].edgeDistanceMap.remove(neighborNP);
						nanoparticles[i].centerDistanceMap.remove(neighborNP);
						nearestNeighbors.remove(neighborNP);
					}
				}
			}*/
			
			
			nanoparticles[i].nearestNeighbors = new Nanoparticle[nearestNeighbors.size()];
			nanoparticles[i].nearestNeighbors = nearestNeighbors.toArray(nanoparticles[i].nearestNeighbors);
			
			nanoparticles[i].nextNearestNeighbors = new Nanoparticle[nextNearestNeighbors.size()];
			nanoparticles[i].nextNearestNeighbors = nextNearestNeighbors.toArray(nanoparticles[i].nextNearestNeighbors);
			
			nanoparticles[i].closeNeighbors = new Nanoparticle[closeNeighbors.size()];
			nanoparticles[i].closeNeighbors = closeNeighbors.toArray(nanoparticles[i].closeNeighbors);
			
			nanoparticles[i].holeCloseNeighbors = new Nanoparticle[closeHoleNeighbors.size()];
			nanoparticles[i].holeCloseNeighbors = closeHoleNeighbors.toArray(nanoparticles[i].holeCloseNeighbors);
			
			nanoparticles[i].neckedNeighbors = new Nanoparticle[neckedNeighbors.size()];
			nanoparticles[i].neckedNeighbors = neckedNeighbors.toArray(nanoparticles[i].neckedNeighbors);
			
			totalnn += nearestNeighbors.size();
			total_neckedNeighbors += neckedNeighbors.size();

		}
		System.out.println("total number of close neighbor is "+totalcn);
		System.out.println("total number of nearest neighbor is "+totalnn);
		System.out.println("total number of close hole neighbors is "+totalhcn);
		System.out.println("total number of necked neighbors is " + total_neckedNeighbors);
		
		average_spacing = total_spacing/total_spaces;
	}
    
    private void throw_electron() {
		int NP, trials=0, maxTrials=100;
		Nanoparticle targetNP;
    	
    	// initialize an electron
    	Electron e = new Electron();
		
    	boolean success = false;
    	
    	while(success != true && trials < maxTrials){
        	// randomly select a NP
    		NP = rng.nextInt(nanoparticles.length);
    		targetNP = nanoparticles[NP];
    		for(int i=0; i<nbnd; i++){
    			if(targetNP.occupationCB[i] < targetNP.occupationMAX[i]){
    				// put electron onto ith orbital on targetNP
    				e.setHost(targetNP, i);
    				targetNP.add_electron(e, i);
    				
    				double[] destination = {targetNP.x, targetNP.z};
    				e.visitedNPs.add(destination);
    				success = true;
    			}
    		}
			trials++;
			
    	}
    	
    	if(!success)
    		System.out.println("Not able to throw electron within max trials!");
	}
    
    private void throw_hole() {
		int NP, trials=0, maxTrials=100;
		Nanoparticle targetNP;
    	
    	// initialize a hole
    	Hole h = new Hole();
		
    	boolean success = false;
    	
    	while(success != true & trials < maxTrials){
        	// randomly select a NP
    		NP = rng.nextInt(nanoparticles.length);
    		targetNP = nanoparticles[NP];
    		for(int i=0; i<nbnd; i++){
    			if(targetNP.occupationVB[i] < targetNP.occupationMAX[i]){
    				// put hole onto ith orbital on targetNP
    				h.setHost(targetNP, i);
    				targetNP.add_hole(h, i);
    				success = true;
    			}
    		}
			trials++;
    	}
    	if(!success)
    		System.out.println("Not able to throw hole within max trials!");
	}
    
    private void throw_exciton() {
    	System.out.println("Threw exciton!");
		int NP, trials=0, maxTrials=100;
		Nanoparticle targetNP;
    	
    	// initialize an electron
    	Electron e = new Electron();
    	Hole h = new Hole();
		
    	boolean success = false;
    	
    	while(success != true && trials < maxTrials){
        	// randomly select a NP
    		NP = rng.nextInt(nanoparticles.length);
    		targetNP = nanoparticles[NP];
    		for(int i=0; i<nbnd; i++){
    			if(targetNP.occupationCB[i] + 1 < targetNP.occupationMAX[i]){
    				// put electron and hole onto ith orbital on targetNP
    				e.setHost(targetNP, i);
    				h.setHost(targetNP, i);
    				targetNP.add_electron(e, i);
    				targetNP.add_hole(h, i);
    				success = true;
    				
    				//also update the events on the nanoparticle
    				targetNP.updateEvents(true, this);
    			}
    		}
			trials++;
			
    	}
    	
    	if(!success)
    		System.out.println("Not able to throw exciton within max trials!");
	}
    
    // for both holes and electron hopping. NOT THOROUGHLY TESTED FOR HOLES
    private void initializeEvents() {
    	rateOnSample = 0.0;
		for(Nanoparticle nanops : nanoparticles){
			for(int band=0; band<nbnd; band++){
				for(Electron electron : nanops.electronsOnNP[band]){
					// look for events and calculate the system total rate
					rateOnSample += electron.lookForEvents(this);

				}
				
				for(Hole hole : nanops.holesOnNP[band]){
					// look for events and calculate the system total rate
					rateOnSample += hole.lookForEvents(this);
					System.out.println("Have holes!");

				}
			}	
		}
	}

    
    // update events, double NP version.
    public void updateEvents(Nanoparticle NP1, Nanoparticle NP2) {
    	// zero current update cycle
    	////currentUpdate.clear();
		// update targetNP and its neighbors
    	NP1.updateEvents(true, this);
    	
    	NP2.updateEvents(true, this);
	}
    
    
    // not thoroughly tested for holes!!
    public HoppingEvent searchHoppingEvent()  {
    	double targetRate = rateOnSample*rng.nextDouble();
    	double currentRate = 0.0 ;
    	
		// loop over all nanoparticles
		for(Nanoparticle nanop : nanoparticles){
			// check each band
			for(int band=0; band<nbnd; band++){
				// loop over electrons sitting on the band
				for(Electron e: nanop.electronsOnNP[band]){
					// finally loop over events associated with the electron
					for(HoppingEvent event : e.hoppings){
						// update latest event to return
						// add current event rate
						currentRate += event.rate;

						if(currentRate>=targetRate){
    						latestEvent = event;
							return latestEvent;
							}
						}
					}
				
				for(Hole h: nanop.holesOnNP[band]){
					// finally loop over events associated with the electron
					for(HoppingEvent event : h.hoppings){
						// update latest event to return
						// add current event rate
						currentRate += event.rate;

						if(currentRate>=targetRate){
    						latestEvent = event;
    						System.out.println("found hole event");
							return latestEvent;
						}
					}
					System.out.println("Have holes!");
				}
			}
    	}

		double totalRate = 0.0;
		for(Nanoparticle nanop : nanoparticles){
			// check each band
			for(int band=0; band<nbnd; band++){
				// loop over electrons sitting on the band
				for(Electron e: nanop.electronsOnNP[band]){
					// finally loop over events associated with the electron
					for(HoppingEvent event : e.hoppings){
						// update latest event to return
						// add current event rate
						totalRate += event.rate;
					}
				}
			}
		}
		System.out.println(totalRate);
		System.out.println("target rate" + targetRate +" less than " + rateOnSample);
    	System.out.println(targetRate<currentRate);
    	if(latestEvent.type=="empty")
    		System.out.println("empty events!"); 
    	return latestEvent;
	}

    
    public double executeEvent(HoppingEvent event, int i) {
    	boolean electronEvent = (event.hostElectron != null);
    	
    	double current=0;
    	
    	if(event.targetNP != previousEvent.sourceNP)
    		event.targetNP.hotness++;
    	
    	if(i >= Configuration.MIN_THRESH_STEPS)
    		//store the electron move in the sourceNP history (where did the electron go)
    		event.sourceNP.add_hop(event.targetNP);
    	
    	if(electronEvent) {
    		Electron eMoving = event.hostElectron;
    		eMoving.move(event.sourceNP, event.sourceOrbital, event.targetNP, event.targetOrbital);
    	}
    	else {
    		Hole hMoving = event.hostHole;
    		hMoving.move(event.sourceNP, event.sourceOrbital, event.targetNP, event.targetOrbital);
    	}
    	
    	// after move
    	if(i >= Configuration.MIN_THRESH_STEPS) {
    		if(electronEvent) {
		    	if(event.targetNP.drain && event.sourceNP.source){
		    		current = -1;
		    		this.sampleCurrent -= 1;
		    	}
		    	
		    	if(event.targetNP.source && event.sourceNP.drain){
		    		current = 1;
		    		this.sampleCurrent +=1 ;
		    	}
    		}
    	}
    	
    	
    	previousEvent = event;
    	return current;
	}
    
    public double illuminationRate() {
    	//under solar illumination by photon flux b_s(E), solar cell absorbs photons at rate (1-R(E))*a(E)*b_s(E)
    	//integrate over energy to obtain total illumination rate, note that only photons above bandgap will be abosrbed
    	double perunitvol = 1.02866e21;
    	return perunitvol*cellx*celly*Math.pow(Constants.ry_to_m,2)*1e3;
    }
    
    public double recombinationRate() {
    	/*Matt Law and Dong Yu Paper Numbers
    	 * capture cross section for e- and e+: 3.8e-18/cm^2
    	 * N_C = N_V = 1e20/cm^3
    	 * E_trap = cbenergy - 80meV
    	 * gapEnergy = .7eV
    	 * trap density N_R = 1e19/cm^3
    	 */
    	
    	//note: we will be assuming the boltzmann approximation holds true
    	double trap_density = 1e-2; //per nm^3
    	double cCS = 3.8e-4; //nm^2
    	
    	//Define gap energy using a random nanoparticle
    	int NP = rng.nextInt(nanoparticles.length);
		Nanoparticle targetNP = nanoparticles[NP];
		double gapEnergy = targetNP.cbenergy[0] - targetNP.vbenergy[0];
		
		//Define trap energy level relative to gap energy
		double trapEnergy = targetNP.cbenergy[0]-.08/Constants.rytoev;
    	
    	//Define Shockley-Reed-Hall Recombination Lifetimes (for holes and electrons)
    	//assume thermal velocities and capture cross sections are the same
    	double thermal_velocity_e = Math.sqrt(2*(3/2*Constants.k_boltzmann_si*temperature*Constants.ry_to_kelvin*300/80)/(8*emass*Constants.electronmass_si))*1e9; // nm/s
    	double thermal_velocity_h = Math.sqrt(2*(3/2*Constants.k_boltzmann_si*temperature*Constants.ry_to_kelvin*300/80)/(8*hmass*Constants.electronmass_si))*1e9; // nm/s
    	double cCS_e = cCS; //elecon capture cross section
    	double cCS_h = cCS; //hole capture cross section
    	double tau_e_SRH = 1/(thermal_velocity_e*cCS_e*trap_density);
    	double tau_p_SRH = 1/(thermal_velocity_h*cCS_h*trap_density);
    	
    	//Define density of states in the conduction and valence bands, per nm^3
    	//double SI_kbT = Constants.k_boltzmann_si*temperature*Constants.ry_to_kelvin;
    	//double h_bar = Constants.h_planck_si/(2*Constants.pi);
    	//double N_V = 2*Math.pow(2*hmass*Constants.electronmass_si*SI_kbT/(2*Constants.pi*Math.pow(h_bar,2)), 1.5)*1e-27; //Jenny Nelson formulation
    	//double N_C = 2*Math.pow(2*emass*Constants.electronmass_si*SI_kbT/(2*Constants.pi*Math.pow(h_bar,2)), 1.5)*1e-27; //Jenny Nelson formulation
    	double N_V = 2*Configuration.h_degeneracy*packingfraction/((4/3.)*Constants.pi*Math.pow(targetNP.diameter*Constants.bohrtonm/2,3)); //QD formulation
    	double N_C = 2*Configuration.e_degeneracy*packingfraction/((4/3.)*Constants.pi*Math.pow(targetNP.diameter*Constants.bohrtonm/2,3)); //QD formulation
    	
    	//Define intrinsic carrier density (thermal equilibrim density of electrons), used to calculate fermi level
    	double n_i = Math.sqrt(N_V*N_C)*Math.exp(-gapEnergy/(2*temperature));
    	
    	//Define free carrier densities
    	double n = nelec/(cellz_nm*cellx*celly*Math.pow(Constants.bohrtonm,2)); //per nm3
    	double p = nhole/(cellz_nm*cellx*celly*Math.pow(Constants.bohrtonm,2)); //per nm3
    	
    	//Define trap carrier densities (fermi energy equal to trap energy)
    	double n_t = N_C*Math.exp((trapEnergy-targetNP.cbenergy[0])/temperature);
    	double p_t = N_V*Math.exp((targetNP.vbenergy[0]-trapEnergy)/temperature);
    	
    	//Define fermi level (using intrinsic density and band weights)
    	double E_i = 0.5*(targetNP.cbenergy[0]+targetNP.vbenergy[0]) - 3/4*temperature*Math.log(emass/hmass);
    	double fermi_energy = temperature*Math.log(n/n_i)+E_i;
    	//double fermi_energy = temperature*Math.log(n/N_C) + targetNP.cbenergy[0];
    	
    	//double n_check = N_C*Math.exp((fermi_energy-targetNP.cbenergy[0])/temperature); //check for self-consistency
    	//double p_check = N_V*Math.exp((-fermi_energy+targetNP.vbenergy[0])/temperature); //check for self-consistency
    	//if n != n_check, or p != p_check, a problem has developed, and its not consistent

    	//electron lifetime (muli-shallow trap SRH)
    	double tau_e_top = tau_e_SRH*(p+p_t) + tau_p_SRH*(n+n_t+trap_density/(1+n/n_t));
    	double tau_e_bottom = n + p + trap_density/((1+n/n_t)*(1+n_t/n)); 
    	double tau_e = tau_e_top/tau_e_bottom;
    	
    	//hole lifetime (multi-shallow trap SRH)
    	double tau_p_top = tau_p_SRH*(n+n_t) + tau_e_SRH*(p+p_t+trap_density/(1+p/p_t));
    	double tau_p_bottom = n + p + trap_density/((1+p/p_t)*(1+p_t/p)); 
    	double tau_p = tau_p_top/tau_p_bottom;
    			
    	//direct recombination rate
    	double baseElectronDensity = 0.0;
    	double tau_p_direct = 1/(temperature*baseElectronDensity);
    	
    	//double holeCC = 1/(1+2*n_i_per_nm3/n_0*Math.cosh((E_trap-fermiLevel)/temperature)); //hole capture coefficient
    	//holeCC = thermal_velocity*cCS*holeCC;
    	
    	//total recombination rate
    	double tau_p_tot = 1/(1/tau_p_direct + 1/tau_p_SRH);
    	
    	//double alt_holedens = Constants.pi/2*Math.pow(4*2*Constants.electronmass_si*.05/Math.pow(Constants.h_planck_si,2), 3/2)*Math.sqrt(targetNP.cbenergy[0]*Constants.rytojoule);
    	
    	System.out.println("Gap Energy is: " + gapEnergy*Constants.rytoev);
    	System.out.println("Thermal velocity is: " + thermal_velocity_e);
    	System.out.println("Density of hole states is: " + N_V*1e21);
    	System.out.println("E_f - E_v: " + (fermi_energy-targetNP.vbenergy[0])*Constants.rytoev);
    	System.out.println("E_c - E_f: "+ (-fermi_energy+targetNP.cbenergy[0])*Constants.rytoev);
    	System.out.println("Electron lifetime: " + tau_e);
    	System.out.println("Hole lifetime: " + tau_p);
    	
    	//return (n-n_0)/tau_p_tot;
    	//have lifetimes, now need recombination rates!
    	//R = (p-p0)/tau_e for electrons
    	//R = (n-n0)/tau_p for holes
    	//n0 and p0 are the thermal equilibrim densities
    	
    	//assume sample is n-type, general recombination rate
    	//double n0 = trap_density;
    	//double p0 = Math.pow(n_i, 2)*trap_density;
    	//double recombination_rate = (p-p0)/tau_p;
    	
    	//generic SRH-style recombination rate formulation
    	double recombination_rate = (n*p - Math.pow(n_i, 2))/(tau_e*(p+p_t) + tau_p*(n+n_t)); //Nelson pg. 108
    	
    	//look into recombination rate at grain boundaries, see Jenny Nelson book pg. 110
    	
    	//double[] lifetimes = {tau_e, tau_p};
    	return recombination_rate;
    }
    
    public void writeElectronDestinations() {
    	int eNum = 0;
    	
    	for(Nanoparticle nanops : nanoparticles){
			for(int band=0; band<nbnd; band++){
				for(Electron electron : nanops.electronsOnNP[band]) {
					String title = "Electron Movements" + "\\" + "Electron" + eNum + ".txt";
					
					PrintWriter writer;
					try {
						writer = new PrintWriter(title);
					    
						for(double[] destination : electron.visitedNPs) {
							double thisx = destination[0]/Constants.nmtobohr;
							double thisz = destination[1]/Constants.nmtobohr;
							String toAdd = thisx + " " + thisz;
							writer.println(toAdd);
						}
			    		writer.close();
					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					eNum += 1;
				}
			}	
		}
	    	
    	
    }
 
    public void writeElectronHoppingVectors() {
    	String title = "Triclinic Lattice Tests" + "\\" + "ElectronVectorField" + ".txt";
		PrintWriter writer;
		try {
			writer = new PrintWriter(title);
	    	for(int j=0; j<nanoparticles.length; j++){
	    		Nanoparticle thisNanoparticle = nanoparticles[j];
	    		double deltaX = 0;
	    		double deltaZ = 0;
	    		Integer totalHops = 0;
	    		
	    		for ( Map.Entry<Nanoparticle, Integer> entry : thisNanoparticle.electronHopMap.entrySet()) {
	    		    Nanoparticle key = entry.getKey();
	    		    Integer value = entry.getValue();
	    		    
	    		    //want the vector for the electron move from this nanoparticle to the key nanoparticle, since electron is hopping in this direction
	    		    //Double[] vector = thisNanoparticle.npnpVector(key, this);
	    		    Double[] vector = thisNanoparticle.npnpVector(key, this);
	    		    
	    		    //normalize vector in x-z so it's a unit vector
	    		    double vectorMag = Math.sqrt(vector[0]*vector[0] + vector[2]*vector[2]);
	    		    
	    		    if(Double.compare(vector[1],0) == 0) {
		    		    deltaX += vector[0]*value/(vectorMag*Constants.nmtobohr);
		    		    deltaZ += vector[2]*value/(vectorMag*Constants.nmtobohr);
		    		    totalHops += value;
	    		    }
	    		}
	    		
	    		//normalize the vector
	    		//double hopDist = Math.sqrt(deltaX*deltaX + deltaZ*deltaZ);
	    		//deltaX = deltaX/totalHops;
	    		//deltaZ = deltaZ/totalHops;
	    		
	    		//get average values
	    		//deltaX = deltaX/totalHops;
	    		//deltaZ = deltaZ/totalHops;
	    		//System.out.println("Nanoparticle: " + j + " had hopping vector: (" + deltaX+", " + deltaZ + ")");
	    		
	    		double thisx = thisNanoparticle.x/Constants.nmtobohr;
	    		double thisy = thisNanoparticle.y/Constants.nmtobohr;
	    		double thisz = thisNanoparticle.z/Constants.nmtobohr;
	    		
	    		String toAdd = thisx + " " + thisy + " " + thisz + " " + totalHops + " " + deltaX + " " + deltaZ;
	    		writer.println(toAdd);
	    	}
    		writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    public void writeElectronAverages() {
    	String currentDirectory = System.getProperty("user.dir");
    	String folder = "New NP necking project" + "\\" + "ElectronAverages" + "\\" + "Disorder" + folderName.substring(folderName.length() - 5) + "\\" + "T" + Math.round(temperature*Constants.ry_to_kelvin*100.0)/100.0;
    	File proposedFolder = new File(currentDirectory +"\\" + folder);
    	System.out.println(currentDirectory + "/" + folder);
    	String title = folder + "\\" + "Sample_" + sample_number+".txt";
    	if(!proposedFolder.exists()) {
    		//System.out.println("Made directory!");
    		boolean madedir = proposedFolder.mkdirs();
    		System.out.println(madedir);
    	}
		PrintWriter writer;
		try {
			writer = new PrintWriter(title);
			writer.println("0 1 2 3");
			int steps = Configuration.STEPS - Configuration.MIN_THRESH_STEPS;
	    	writer.println(totalZeroElectronNanoparticles/steps + " " + totalOneElectronNanoparticles/steps + " "+ totalTwoElectronNanoparticles/steps + " " + totalThreeElectronNanoparticles/steps);
    		writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    public double simulation() throws FileNotFoundException{
    	
    	
    	
    	long l;
    	Nanoparticle source, target;
    	
    	l = System.nanoTime();
    	
        initializeEvents();
        double fractionalIncidentPhotons = 0.0;
        double startTime = 0.0;
        double startCurrent = 0.0;
        double startSteps = Configuration.MIN_THRESH_STEPS;
        double iteration = 1;
        System.out.println("rateOnSample is: " + rateOnSample);
        
//        String title = "TriDENS Results" + "\\" + "MobilityLogSample" + sample_number + ".txt";
//		PrintWriter writer = new PrintWriter(title);
//		writer.println("Step Mobility(cm^2/Vs)");
    	
    	for(int i=0; i<Configuration.STEPS; i++){
    		
    		//System.out.println();
    		//System.out.println("step "+i);
            currentEvent = new HoppingEvent();
    		currentEvent = searchHoppingEvent();

    		source = currentEvent.sourceNP;
    		target = currentEvent.targetNP;
    		
    		executeEvent(currentEvent, i);
    		double ry_timestep = -Math.log(rng.nextDouble())/rateOnSample; //time-step in Ry units
    		if(i >= Configuration.MIN_THRESH_STEPS) {
    			elapsedTime += ry_timestep*Constants.ry_ps; //elapsed time in ps
    			for(int j=0; j<nanoparticles.length; j++){
    	    		int nElectrons = nanoparticles[j].occupationTotalElectron;
    	    		if(nElectrons == 0) totalZeroElectronNanoparticles += 1.0;
    	    		if(nElectrons == 1) totalOneElectronNanoparticles += 1.0;
    	    		if(nElectrons == 2) totalTwoElectronNanoparticles += 1.0;
    	    		if(nElectrons == 3) totalThreeElectronNanoparticles += 1.0;
    			}
    		}
    		
//    		if(i - startSteps > 100000) {
//    			double thisTime = elapsedTime - startTime;
//    			double thisCurrent = sampleCurrent - startCurrent;
//    			double thismobility = thisCurrent*cellz_nm*cellz_nm*.01/(nelec*voltage*Constants.ry_volt*thisTime); //cm^2/Vs
//    			//System.out.println("current mobility is: " + thismobility + " steps: " + (i - startSteps) + " current: " + (sampleCurrent-startCurrent));
//    			//System.out.println("current mobility is: " + thismobility);
//    			writer.println(iteration + " " + thismobility);
//    			startSteps = i;
//    			startCurrent = sampleCurrent;
//    			startTime = elapsedTime;
//    			iteration++;
//    		}
    		
    		//initializeEvents();
    		updateEvents(source, target);
//    		double totalRate = 0.0;
//    		for(Nanoparticle nanop : nanoparticles){
//    			// check each band
//    			for(int band=0; band<nbnd; band++){
//    				// loop over electrons sitting on the band
//    				for(Electron e: nanop.electronsOnNP[band]){
//    					// finally loop over events associated with the electron
//    					for(HoppingEvent event : e.hoppings){
//    						// update latest event to return
//    						// add current event rate
//    						totalRate += event.rate;
//    					}
//    				}
//    			}
//    		}
//    		System.out.println("total rate is: " + totalRate);
//    		System.out.println("Rate on sample is: " + rateOnSample);
    		
    		//make map of electron flow by storing all hops
    		
    		//System.out.println(recombinationRate());
    		//double[] recomRate = recombinationRate();
    		//double e_rate = recomRate[0];
    		//double h_rate = recomRate[1];
    		
    		//do we want to include recombination and illumination in the kinetic monte carlo rate table?
    		
    		//System.out.println("Illumination rate: " + illuminationRate());
    		/*
    		double SI_timestep = ry_timestep*Constants.ry_ps*1.0e-12;
    		fractionalIncidentPhotons += SI_timestep*illuminationRate();
    		
    		int numberIncidentPhotons = (int)fractionalIncidentPhotons;
    		//throw excitons depending on illumination rate, the throw_exciton method will update events as electrons are added
    		for(int iterator = 0; iterator < numberIncidentPhotons; iterator++) {
    			throw_exciton();
    			System.out.println("number photons: " + numberIncidentPhotons);
    			System.out.println("Throwing exciton");
    		}
    		
    		if(fractionalIncidentPhotons > 1.0) {
    			System.out.println("Reset number of photons");
    			fractionalIncidentPhotons = 0.0;
    		}
    		*/
    		//System.out.println("SI time: " + SI_time);
    		
    		//Trying to implement grand canonical ensemble
    		/*if(i%20000==0){
    			
    			System.out.println(i);
    			source.add_electron_try(this);
    
    		}*/
    		
    		
    		
    		
    		/*/debug section
    		System.out.println(currentEvent+" "+currentEvent.sourceNP+" "+currentEvent.targetNP);
    		numberEvents = 0;
    		for(Nanoparticle nanop : nanoparticles){
    			for(Electron e:nanop.electronsOnNP[0]){
    				numberEvents += e.hoppings.size();
    			}
    		}
    		System.out.println("number of events "+ numberEvents);
    		
    		if(i%20000==0){
    			
    			System.out.println(i);
    			//System.out.println(currentEvent.sourceNP);
    			System.out.println(sampleCurrent);
    		}
    		*/
    		
    		
    		
    		// end of debug section
    		

    	}
    	//writer.close();
		System.out.println("sample current is: " + sampleCurrent);
		System.out.println("elapsed time is : " + elapsedTime );
		
		//writeElectronHoppingVectors();
		//writeElectronDestinations();
		double mobility = sampleCurrent*cellz_nm*cellz_nm*.01/(nelec*voltage*Constants.ry_volt*elapsedTime); //cm^2/Vs
		System.out.println("mobility is: " + mobility);

    	l = System.nanoTime() - l;
        System.out.println("iteration took " + l/1000000000 + "s");
	        
        writeElectronAverages();
        //return mobility;
        return mobility;
    }
 
    
    

	public static void main(String[] args) throws FileNotFoundException {

		// Necessary parameters
		//(INT_t steps, INT_t sample_no, INT_t e_number, INT_t h_number, INT_t nanops_number, 
		//FLOAT_t large_ratio, FLOAT_t temp, FLOAT_t thr, str features, FLOAT_t voltage_ratio, output):
        Map<String, Object> params = new HashMap<>();
        
        params.put("nelec", 100);
    	params.put("nholes", 0);
    	params.put("expected_nnanops", 400);
    	params.put("sample_no", 0);
    	params.put("feature", "mobility");
    	params.put("temperature", 80.0);
    	params.put("thr", 0.0);
    	params.put("np_diam", 3.0);
        
        Sample newsample = new Sample(params);
        
        newsample.simulation();
        
        
        
        
        //System.out.println(newsample.getElectronMass());
        //System.out.println(newsample.getFWHM()*Constants.rytoev);
        //Nanoparticle[] nanoparticles;
        //nanoparticles = setupNanoparticles(newsample);
        
        //System.out.println(nanoparticles[0].x/Constants.nmtobohr);
        //System.out.println(nanoparticles[399].getCB1());
	}
}
    

