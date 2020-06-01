package util;

import com.sun.org.apache.bcel.internal.generic.NEW;

public class Configuration {
	
	// System
	public static final int STEPS = 250000; //total number of steps
	public static final int MIN_THRESH_STEPS = 100000; //number of steps before values are recorded
	public static final double ligandDC = 2.0; // ligand dc. EDT: 2.0. Alumina: 9.0. Used 9.0 for Moule
	public static final double npDC = 22.0; //bulk dielectric constant. 22.9 more accurate. If T dependent, bulk@300K.
	public static final double ligandLength = 0.50 ; // in nm, was 0.2 for Moule project //was 0.45 for tridens
	public static final double emass = 0.05/2.0; //effective mass for electrons 
	public static final double hmass = 0.05/2.0; //effective mass for holes 
	public static final String hoppingMechanism = "ma"; // 'ma' : miller-abrahams, 'marcus' : marcus
	public static final String capacitanceModel = "delerue"; //"delerue" or "zunger"
	public static final String effectiveMedium = "shklovskii local"; //"mg" or "shklovskii global" or "shklovskii local"
	public static final String poissonSolver = "none";
	public static final int[] rngSeed = new int[]{1,2,3,4};
	
	// Lattice
	public static final boolean PERX = true; //whether the structure is periodic in hopping in x
	public static final boolean PERY = true; //whether the structure is periodic in hopping in y
	public static final boolean PERZ = true; //whether the structure is periodic in hopping in z
	public static final double latticeAngleGamma = 90.0*Math.PI/180.0; //lattice basis angle in radians //96 for Moule //99 for tridens
	public static final double latticeAngleAlpha = 90.0*Math.PI/180.0; //lattice basis angle in radians //103 for Moule //99 for tridens
	public static final double latticeAngleBeta = 90.0*Math.PI/180.0; //lattice basis angle in radians //95 for Moule //99 for tridens
	public static final String latticeStructure = "cubic"; //"triclinic" or "cubic" or "Moule". Was triclinic for Moule project
	public static final double latticeBasisLength = 7.5; //in nm
	public static final String transportDirection = "z"; //"x", "y", or "z"
	
	//Necking NOT USED CURRENTLY
	public static final double neckRadiusMean = 3.0; //in nm
	public static final double neckRadiusSigma = 0.5; //in nm
	public static final double neckOccurenceMean = 3.0; //number of necked neighbors
	public static final double neckOccurenceSigma = 0.5;
	public static final double neckDistThr = 1.0; //in nm, the edge-to-edge distance cutoff for necking to occur
	
	public static final double edgeDistThr = 3.5; //in nm
		
	// NP
	public static final double near_n_distThr = 3.0; //in nm, for cubic
	public static final double nextn_n_distThr = 5.0; //in nm, for cubic
	//public static final double distThr = 4.5; //in nm, for triclinic
	//public static final double distThr = 1.27;
	//public static final double distThr = 3.2;
	public static final double bradii = 47.0;
	public static final double reorgenergy = 0.05;
	public static final double jumpFreq = 1.0/7.8; //nu, tune to adjust mobility values. Was 0.7 for comparison to Kang, 1.0 normal
	public static final double marcusprefac2 = 1.0;
	public static final double capacitance = 0.05; //capacitance is used if capacitance_model = 'zunger'
	public static final double screeningFactor = 1.0; //screening factor to change charging energy
	public static final int nBands = 1; //should add hole and electrons bands separately, 2 e bands, 1 h band
	public static final double metallicPrefactor = 1.0; //metallic jump frequency. Should change to be T and material dependent
	public static final String diameter = "3.0nm/";
	public static final int e_degeneracy = 8;
	public static final int h_degeneracy = 8;
	public static final double epsilon_per_K = -.01; //change in dielectric constant per K
	
	// Boolean control
	public static final boolean twoLayer = false;
	public static final boolean biModal = false;
	public static final boolean pennModel = false;
	public static final boolean lcapacitance0 = true; //whether to add polarization term
	
	
}
