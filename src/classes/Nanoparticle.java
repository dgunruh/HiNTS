package classes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

import util.Configuration;
import util.Constants;

public class Nanoparticle {
	
	public double diameter, radius, ligandlength, x, y, z;
	public int id;
	
	int nn; //number of nearest neighbors
	Nanoparticle[] nearestNeighbors;  //nearest neighbor array
	
	int nnn; //number of next nearest neighbors
	Nanoparticle[] nextNearestNeighbors;
	
	int cn; //number of closest neighbors
	Nanoparticle[] closeNeighbors;  //close neighbor array
	
	int hcn;
	Nanoparticle[] holeCloseNeighbors;  //close neighbor array
	
	int neckedn; //number of necked neighbors
	Nanoparticle[] neckedNeighbors; //necked nanoparticle array
	
	int nelectrons; // number of electrons on board
	LinkedList<Electron>[] electronsOnNP;   //use LinkedList here because we will do a lot of electron deletion and insertion.
	
	int nholes; // number of holes on board
	LinkedList<Hole>[] holesOnNP;
	
	
	HashMap<Nanoparticle, Double> edgeDistanceMap; //distance map to neighbors, edge to edge
	HashMap<Nanoparticle, Double> centerDistanceMap; //distance map to neighbors, center to center
	HashMap<Nanoparticle, Integer> electronHopMap; //number of electron hops from each nanoparticle
	HashMap<Nanoparticle, Double> neckRadiusMap; //radius map of necked nanoparticles

	
	long hotness;
	
	double[] cbenergy = new double[2];
	double cbenergy1, cbenergy2;
	
	double[] vbenergy = new double[2];
	double vbenergy1, vbenergy2;
	
	double dcin, dcout;
	
	double selfenergy0, selfenergy;
	double hselfenergy0, hselfenergy;
	
	public boolean inCluster, source, drain;
	int clusterIndex; 
	
	int[] occupationCB, occupationVB, occupationMAX;
	int occupationTotalElectron, occupationTotalHoles;
	//int[] electronIndex, holeIndex;
	
	public Nanoparticle(double xcoord, double ycoord, double zcoord, double diameter_input, int id_value, Sample sample) {
		// TODO Auto-generated constructor stub
		x = xcoord*Constants.nmtobohr;
		y = ycoord*Constants.nmtobohr;
		z = zcoord*Constants.nmtobohr;
		double x0 = Configuration.latticeBasisLength*Constants.nmtobohr/2.0;
		double y0 = Configuration.latticeBasisLength*Constants.nmtobohr/2.0;
		double z0 = Configuration.latticeBasisLength*Constants.nmtobohr/2.0;
		diameter = diameter_input*Constants.nmtobohr;// - 2.0*sample.ligandlength; //question about this
		radius = diameter/2;
		id = id_value;
		set_cbenergy();
		set_vbenergy();
		set_occupation(sample);
		
		//Define whether it is a source or a drain
		switch(Configuration.transportDirection) {
		case "z":
			switch(sample.latticeStructure){
			case "triclinic":
				if(!sample.reverse_the_sample) {
					//The sample does not need to be reversed. It has either been flipped or is in its standard configuration
					//calculate components of lattice vector a_3
					double x_3 = Configuration.latticeBasisLength*(Math.cos(sample.latticeAngleBeta)-Math.cos(sample.latticeAngleAlpha)*Math.cos(sample.latticeAngleGamma))/Math.sin(sample.latticeAngleGamma);
					double z_3,m_edge;

					if(sample.sample_number%2 == 0) {
						z_3 = Configuration.latticeBasisLength*Math.cos(sample.latticeAngleAlpha);
						m_edge = 1/Math.tan(sample.latticeAngleGamma); //equation is z-z0 = m*(x-x0
					}
					else {
						z_3 = -Configuration.latticeBasisLength*Math.cos(sample.latticeAngleAlpha);
						m_edge = -1/Math.tan(sample.latticeAngleGamma); //equation is z-z0 = m*(x-x0
					}


					double y_3 = Math.sqrt(Math.pow(Configuration.latticeBasisLength,2.0)-Math.pow(x_3,2.0)-Math.pow(z_3,2.0));
					int y_layer = (int) Math.round((y-y0)/(y_3*Constants.nmtobohr));
					//calculate how the x and z components of the origin of that layer are shifted with respect to the first layer
					double z_mod = y_layer*z_3*Constants.nmtobohr;
					double x_mod = y_layer*x_3*Constants.nmtobohr;


					//calculate the left and right electrode locations for the current NP
					double z_left = m_edge*(x-x_mod - x0) + z_mod + z0 - Configuration.ligandLength*Constants.nmtobohr/2.0;
					double z_right = sample.cellz + m_edge*(x-x_mod - x0) + z_mod - z0 + Configuration.ligandLength*Constants.nmtobohr/2.0;

					//if the edge of the NP overlaps with the electrode, it is considered to be a source/drain
					if(z-radius <= (z_left + sample.near_n_dist_thr)){
						source = true;
						sample.sources.add(this);
					}
					if(z+radius >= (z_right - sample.near_n_dist_thr)){
						drain = true;
						sample.drains.add(this);
					}
				}
				else {
					//Sample is to be flipped internally, so sources and drains are reversed
					if(z-radius <= sample.near_n_dist_thr){
						drain = true;
						sample.drains.add(this);
					}
					if(z+radius >= (sample.cellz - sample.near_n_dist_thr)){
						source = true;
						sample.sources.add(this);
					}
				}
				break;
			case "cubic":
				if(z-radius <= sample.near_n_dist_thr){
					source = true;
					sample.sources.add(this);
				}
				if(sample.cellz - (z+radius) <= sample.near_n_dist_thr){
					drain = true;
					sample.drains.add(this);
				}
				break;

			case "Moule":
				double z_left_edge = sample.leftmost_NP_z;
				double z_right_edge = sample.rightmost_NP_z;
				if(z-radius <= (z_left_edge + sample.near_n_dist_thr)){
					source = true;
					sample.sources.add(this);
				}
				if(z+radius >= (z_right_edge - sample.near_n_dist_thr)){
					drain = true;
					sample.drains.add(this);
				}
				break;
		}
		break;
		case "x":
			if(!sample.reverse_the_sample) {
				if(x-radius <= sample.near_n_dist_thr){
					source = true;
					sample.sources.add(this);
				}
				if(sample.cellx - (x+radius) <= sample.near_n_dist_thr){
					drain = true;
					sample.drains.add(this);
				}
			}
			else {
				if(x-radius <= sample.near_n_dist_thr){
					drain = true;
					sample.drains.add(this);
				}
				if(sample.cellx - (x+radius) <= sample.near_n_dist_thr){
					source = true;
					sample.sources.add(this);
				}
			}
			break;

		case"y":
			if(!sample.reverse_the_sample) {
				if(y-radius <= sample.near_n_dist_thr){
					source = true;
					sample.sources.add(this);
				}
				if(sample.celly - (y+radius) <= sample.near_n_dist_thr){
					drain = true;
					sample.drains.add(this);
				}
			}
			else {
				if(y-radius <= sample.near_n_dist_thr){
					drain = true;
					sample.drains.add(this);
				}
				if(sample.celly - (y+radius) <= sample.near_n_dist_thr){
					source = true;
					sample.sources.add(this);
				}
			}
			break;
		}
		
		electronsOnNP = new LinkedList[2];
		holesOnNP = new LinkedList[2];
		for(int i=0; i<2; i++){
			electronsOnNP[i] = new LinkedList<Electron>();
			holesOnNP[i] = new LinkedList<Hole>();
		}
		
	}
	
	
	private void set_cbenergy(){
		 // Kand and Wise, d in nm, energy in eV
		 // cbenergy=evtory*(0.238383 +4.28063 *diameter**(-1.87117 ))
		 // cbenergy2=0.254004 +8.0018 *diameter**(-1.82088 )
		 // localdiam=diameter*bohrtonm
         cbenergy1 = Constants.evtory*(-4.238383 + 4.28063*Math.pow((diameter*Constants.bohrtonm), -1.87117));
		 cbenergy[0] = Constants.evtory*(-4.238383 + 4.28063*Math.pow((diameter*Constants.bohrtonm), -1.87117));
		 cbenergy2 = Constants.evtory*(-4.254004 + 8.0018*Math.pow((diameter*Constants.bohrtonm), -1.82088));
	     cbenergy[1] = Constants.evtory*(-4.254004 + 8.0018*Math.pow((diameter*Constants.bohrtonm), -1.82088));
	}
	
	private void set_vbenergy(){
	     // Kand and Wise, d in nm, energy in eV
	     // localdiam=diameter*bohrtonm
	     // vbenergy=evtory*(-0.220257 -5.23931 *diameter**(-1.85753 ))
		 vbenergy1 = Constants.evtory*(-4.728265 - 5.23931*Math.pow((diameter*Constants.bohrtonm), -1.85753));
		 vbenergy[0] = Constants.evtory*(-4.728265 - 5.23931*Math.pow((diameter*Constants.bohrtonm), -1.85753));
	}
	
	private void set_occupation(Sample sample){
		occupationCB = new int[sample.nbnd];
		occupationVB = new int[sample.nbnd];
		occupationMAX = new int[sample.nbnd];
		for(int i=0; i<sample.nbnd; i++){
			occupationMAX[i] = sample.e_degeneracy;  // for symmetric CB and VB only
		}
		occupationTotalElectron = 0;
		occupationTotalHoles = 0;
	}
	
	
	/**
	@param otherNP, sample, if periodic in z, if center to center
	@return NP-NP distance
	@throws what kind of exception does this method throw
	*/
	public double npnpdistance(Nanoparticle otherNP, Sample sample, boolean connected_z, boolean CenterToCenter) {
		
		double deltax = 0.0, deltay,deltaz = 0.0, distSquared;
		double cellx = sample.cellx;
		double celly = sample.celly;
		double cellz = sample.cellz;
		if (Configuration.PERX) {
			deltax = Math.min(Math.min(Math.abs(this.x-otherNP.x), Math.abs(this.x-(otherNP.x-cellx))), Math.abs(this.x-(otherNP.x+cellx)));
			if (Configuration.latticeStructure == "triclinic") {
				if (Math.abs(deltax - Math.abs(this.x - otherNP.x)) > .1) {
					double z_2 = Configuration.latticeBasisLength*Math.cos(sample.latticeAngleGamma);
					int y_layer_other = (int) otherNP.id/(Configuration.x_nps*Configuration.z_nps);
					int y_layer_this = (int) this.id/(Configuration.x_nps*Configuration.z_nps);
					
					int thisx_int = (int) (this.id - (y_layer_this*Configuration.x_nps*Configuration.z_nps))/Configuration.z_nps;
					int otherx_int = (int) (otherNP.id - (y_layer_other*Configuration.x_nps*Configuration.z_nps))/Configuration.z_nps;
					
					//First correct for the distance that the deltaz calc will give
					deltaz -= Math.abs(otherx_int - thisx_int)*z_2*Constants.nmtobohr;
					
					int delta_NP_xlayers = Math.abs(Configuration.x_nps - Math.abs(otherx_int - thisx_int));
					deltaz += delta_NP_xlayers*z_2*Constants.nmtobohr;
				}
			}
		}
		else {
			deltax = this.x - otherNP.x;
		}
		
		if(Configuration.PERY)
			deltay = Math.min(Math.min(Math.abs(this.y-otherNP.y), Math.abs(this.y-(otherNP.y-celly))), Math.abs(this.y-(otherNP.y+celly)));
		else
			deltay = this.y - otherNP.y ;
		
		if(connected_z)
			deltaz += Math.min(Math.min(Math.abs(this.z-otherNP.z), Math.abs(this.z-(otherNP.z-cellz))), Math.abs(this.z-(otherNP.z+cellz)));
		else
			deltaz += this.z - otherNP.z ;

		distSquared = Math.pow(deltax, 2.0) + Math.pow(deltay, 2.0) + Math.pow(deltaz, 2.0);

		if(distSquared==0){
			//System.out.println("Same NP!");
			return 0.0;
		}else{
			if(CenterToCenter)
				return Math.pow(distSquared, 0.5);
			else
				return Math.pow(distSquared, 0.5)-(this.radius + otherNP.radius);
		}
		
	}
	
	/**
	@param point, sample
	@return NP-point distance
	@throws no exceptions
	*/
	public double nppointdistance(double[] point, Sample sample) {
		
		double deltax, deltay,deltaz = 0.0, distSquared;
		double cellx = sample.cellx;
		double celly = sample.celly;
		double cellz = sample.cellz;
		
		if(Configuration.PERX) {
			deltax = Math.min(Math.min(Math.abs(this.x-point[0]), Math.abs(this.x-(point[0]-cellx))), Math.abs(this.x-(point[0]+cellx)));
			if (Configuration.latticeStructure == "triclinic") {
				if (Math.abs(deltax - Math.abs(this.x - point[0])) > .1) {
					double x_2 = Configuration.latticeBasisLength*Math.sin(sample.latticeAngleGamma);
					double z_2 = Configuration.latticeBasisLength*Math.cos(sample.latticeAngleGamma);
					double x_3 = Configuration.latticeBasisLength*(Math.cos(sample.latticeAngleBeta)-Math.cos(sample.latticeAngleAlpha)*Math.cos(sample.latticeAngleGamma))/Math.sin(sample.latticeAngleGamma);
					double z_3 = Configuration.latticeBasisLength*Math.cos(sample.latticeAngleAlpha);
					double y_3 = Math.sqrt(Math.pow(Configuration.latticeBasisLength,2.0)-Math.pow(x_3,2.0)-Math.pow(z_3,2.0));
					int y_layer_point = (int) Math.round((point[1]-Configuration.latticeBasisLength*Constants.nmtobohr/2.0)/(y_3*Constants.nmtobohr));
					int y_layer_np = (int) this.id/(Configuration.x_nps*Configuration.z_nps);
					
					double x_mod = y_layer_point*x_3;
					
					int npx_int = (int) (this.id - (y_layer_np*Configuration.x_nps*Configuration.z_nps))/Configuration.z_nps;
					double point_deltax = point[0]*Constants.bohrtonm - x_mod;
					
					int pointx_int = (int) Math.round(point_deltax/x_2) - 1;
					
					//System.out.println("Point y layer: " + y_layer_point + " point y: " +point[1]*Constants.bohrtonm+" Point x layer: " + pointx_int + " point x: "+ point[0]*Constants.bohrtonm +"  point deltax: " + point_deltax + " point y: " + point[1]*Constants.bohrtonm);
					//First correct for the distance that the deltaz calc will give
					deltaz -= Math.abs(pointx_int - npx_int)*z_2*Constants.nmtobohr;
					
					int delta_NP_xlayers = Math.abs(Configuration.x_nps - Math.abs(pointx_int - npx_int));
					deltaz += delta_NP_xlayers*z_2*Constants.nmtobohr;
				}
			}
		}
		else
			deltax = this.x - point[0] ;
		
		if(Configuration.PERY)
			deltay = Math.min(Math.min(Math.abs(this.y-point[1]), Math.abs(this.y-(point[1]-celly))), Math.abs(this.y-(point[1]+celly)));
		else
			deltay = this.y - point[1] ;
		
		if(Configuration.PERZ)
			deltaz += Math.min(Math.min(Math.abs(this.z-point[2]), Math.abs(this.z-(point[2]-cellz))), Math.abs(this.z-(point[2]+cellz)));
		else
			deltaz += this.z - point[2] ;
			
		distSquared = Math.pow(deltax, 2.0) + Math.pow(deltay, 2.0) + Math.pow(deltaz, 2.0);

		return Math.pow(distSquared, 0.5);
	}
		
	
	/**
	@param otherNP, sample, if periodic in z, if center to center
	@return NP-NP distance
	@throws what kind of exception does this method throw
	*/
	//vector to the other nanoparticle
	public Double[] npnpVector(Nanoparticle otherNP, Sample sample) {
		
		double deltax, deltay,deltaz = 0.0;
		double cellx = sample.cellx;
		double celly = sample.celly;
		double cellz = sample.cellz;
		
		if(Configuration.PERX)
			deltax = Math.min(Math.min(Math.abs(this.x-otherNP.x), Math.abs(this.x-(otherNP.x-cellx))), Math.abs(this.x-(otherNP.x+cellx)));
		else
			deltax = otherNP.x - this.x ;
		
		if(Configuration.PERY)
			deltay = Math.min(Math.min(Math.abs(this.y-otherNP.y), Math.abs(this.y-(otherNP.y-celly))), Math.abs(this.y-(otherNP.y+celly)));
		else
			deltay = otherNP.y - this.y ;
		
		if(Configuration.PERZ)
			if(Math.abs(this.z-otherNP.z) <= Math.abs(this.z-(otherNP.z-cellz)))
				deltaz = otherNP.z - this.z;
			
			else 
				deltaz = (otherNP.z - cellz) - this.z;
			
			
			if(Math.abs(deltaz) >= Math.abs(this.z - (otherNP.z+cellz))) 
				deltaz = (otherNP.z + cellz) - this.z;
			
		else
			deltaz = otherNP.z - this.z ;
		
		Double[] vector = {deltax, deltay, deltaz};
		return vector;
	}
	
	private void set_shiftCB(double energyShift){
		cbenergy[0] += energyShift;
	}
	
	private void set_shiftVB (double energyshift){
		vbenergy[0] += energyshift;
	}
	
	public double getCB1() {
		return cbenergy[0];
	}
	
	public int[] getCBoccupation() {
		return occupationCB;
	}
	
	public int[] getMAXocc() {
		return occupationMAX;
	}
	
	// simply add one electron
	public void add_electron(Electron e, int band){
		//System.out.println("adding to "+this);
		//e.hostNP = this;
		electronsOnNP[band].add(e);
		occupationCB[band]++;
		occupationTotalElectron++;
		
	}
	
	// TODO : implement these methods in the Nanoparticle class
	// used for the grand canonical ensemble
	// works for ONE band only
	public boolean add_electron_try(Sample sample) {
		
		double energy_before, energy_after, energy_diff = 0;
		double probability, debug;
		
		for(int band=0; band<sample.nbnd; band++){
			if(occupationCB[band]<occupationMAX[band]){
				

				//energy before insertion, using only one band
				energy_before = totalEnergy(occupationCB[0], occupationVB[0], 0, sample);
				  
				//energy after insertion
				energy_after = totalEnergy(occupationCB[0]+1, occupationVB[0], 0, sample);
				
				double energy_differ = energy_before - energy_after;
				System.out.println("Energy before - after = " + energy_differ);
				System.out.println("Sample temperature in Ry: " + sample.temperature);
				
				//sample.chemicalPotential = -((cbenergy[0] + 4.4*Constants.evtory) + (cbenergy[0] + vbenergy[0])/2)/sample.nelec + sample.temperature*Math.log(sample.nelec/sample.V_prime);
				//sample.chemicalPotential = +(4.4*Constants.evtory)/sample.nnanops + sample.temperature*Math.log(sample.nelec/sample.V_prime);
				//sample.chemicalPotential = sample.temperature*Math.log(sample.nelec/sample.V_prime);
				//(cbenergy[0] + 4.4*Constants.evtory) + 
				//sample.chemicalPotential = -2*0.14*Constants.evtory;
				//sample.chemicalPotential = (cbenergy[0] + vbenergy[0])/2;
				System.out.println("Sample chemical potential is: " + sample.chemicalPotential + " in Ry. It's: " + sample.chemicalPotential*Constants.rytoev +" in eV.");
				double gapEnergy = cbenergy[0] - vbenergy[0];
				double intrinsic_density_per_m3 = 2*Math.pow(Constants.k_boltzmann_si*sample.temperature*Constants.ry_to_kelvin/(Math.pow(Constants.h_planck_si,2)/(2*Constants.pi)), 1.5)*Math.pow(4*Configuration.emass*Configuration.hmass*Math.pow(Constants.electronmass_si,2), .75)*Math.exp(-gapEnergy*Constants.rytojoule/(2*sample.temperature*Constants.k_boltzmann_si*Constants.ry_to_kelvin));
				double intrinsic_density_per_nm3 = intrinsic_density_per_m3*10e-27;
				double density_hole_states = 2*Configuration.h_degeneracy*sample.packingfraction/((4/3.)*Constants.pi*Math.pow(diameter*Constants.bohrtonm/2,3));
				//intrinsic_density_per_nm3 = intrinsic_density_per_nm3*10e-10;
				double fermiLevel = sample.temperature*Math.log(density_hole_states/intrinsic_density_per_nm3) + vbenergy[0];
				double lambda = Math.pow(2*Math.PI/(Configuration.emass*sample.temperature), 0.5)*Constants.bohr_radius_si;
				double rel_fermiLevel = fermiLevel - sample.temperature*Math.log(Math.pow(lambda, 3));
				double chemicalPotential = sample.chemicalPotential;
				System.out.println("Sample intrinsic carrier density is: " + intrinsic_density_per_nm3);
				System.out.println("Density of holes states calculated is: " + density_hole_states);
				System.out.println("Fermi Level using this is: " + fermiLevel + " in Rydbergs");
				System.out.println("Relative Fermi Level: " + rel_fermiLevel);
				System.out.println("Fermi Level using this is: " + fermiLevel*Constants.rytoev + " in electron volts");
				System.out.println("Valence band is at: " + vbenergy[0]*Constants.rytoev + " in eV. Conduction band is at: " + cbenergy[0]*Constants.rytoev + " in eV.");
				
				// energy difference
				//energy_diff = rel_fermiLevel + energy_before - energy_after ;  //sample.chemicalPotential
				energy_diff =  energy_before - energy_after;
				
				
				// TODO: working here now 12_30
				// The numbers are HUGE, does not make any sense, need to correct!!!!!!
				probability = Math.min(1, (sample.V_prime/ (sample.nelec+1))*Math.exp(energy_diff/(sample.temperature)));
				//Math.exp(energy_diff/(Constants.k_boltzmann_ry*sample.temperature))
				
				//debug = Math.exp(energy_diff/Constants.k_boltzmann_ry*sample.temperature);
				//debug = energy_diff/(sample.temperature);
				debug = (sample.V_prime/ (sample.nelec+1))*Math.exp(energy_diff/(sample.temperature));
				
				//System.out.println("adding probability is "+probability);
				System.out.println("Debug is " + debug);
				return true;
			}
		}
		
		
		return false;
		
	}
	
	public double calculateOnSiteElectronCharging(int occupationNumber){
		return occupationNumber * (selfenergy0 + (occupationNumber - 1)*selfenergy/2);
	}
	
	public double calculateOnSiteHoleCharging(int occupationNumber) {
		return occupationNumber * (hselfenergy0 + (occupationNumber -1)*hselfenergy/2);
	}
	
	public double calculateOnSiteExcitonCharging(int eOccupation, int hOccupation, Sample sample) {
		//this uses the binding energy for excitons constrained to NPs, calculated under the effective mass approximation
		//See Delerue pg. 111-112
		double C_coul = 1.786; //for approximately spherical nanoparticles
		double U_direct = C_coul*2*Constants.e2/(dcin*diameter); //direct electron-hole interaction
		double U_image = 2*Constants.e2*(dcin - dcout)/(diameter*dcin*dcout);
		double U_exciton = U_direct + U_image;//note: this will be made negative if electrons, and kept positive if holes
		
		return eOccupation * hOccupation * U_exciton; //each majority carrier will have interactions with each minority carrier
	}
	
	// Works for ONE band only !!!!
	private double totalEnergy(int numElectrons, int numHoles, int band, Sample sample) {
		
		double energy = 0;
		
		// total kinetic energy
		//energy += (numElectrons * cbenergy[band] + numHoles * vbenergy[band]);
				
		// external 'horizontal' potential
		//energy += -1.0* Constants.sqrt2 * sample.voltage * this.z * numElectrons ;   
		energy += -1.0* Constants.sqrt2 * (sample.voltage/sample.cellz)*this.z * numElectrons;
		//energy += 1.0* Constants.sqrt2 * sample.voltage * this.z * numHoles ;   

		// onsite charging
		if(numElectrons > 0){
			/* old method
			// add first electron
			energy += (Configuration.lcapacitance0)? selfenergy0 : 0;
			//System.out.println(selfenergy0+" self "+selfenergy);
			// add subsequent electron
			energy += (numElectrons>1) ? (numElectrons)*selfenergy : 0;
			*/
			
			//new method
			energy += calculateOnSiteElectronCharging(numElectrons);
		} 
		if(numHoles > 0){
			/* old method
			// add first electron
			energy += (Configuration.lcapacitance0)? hselfenergy0 : 0;
			// add subsequent electron
			energy += (numHoles>1) ? (numHoles)*hselfenergy : 0;
			*/
			
			//new method
			energy += calculateOnSiteHoleCharging(numHoles);
		}
		
		// TODO: ignoring excitonic effect for the moment, as there should be ZERO holes.
		System.out.println("energy is "+energy);
		return energy;
	}
	
	// TODO
	// used for the grand canonical ensemble
	// try to add one electron to the NP. 
	public boolean remove_electron_try(Sample sample) {
		return false;
	}
	
	// simply remove one electron
	public void remove_electron(Electron e, int band){
		//System.out.println("removing from "+e.hostNP+" "+this);
		electronsOnNP[band].remove(e);
		occupationCB[band]--;
		occupationTotalElectron--;
	}
	
	
	public void add_hole(Hole h, int band){
		holesOnNP[band].add(h);
		occupationVB[band]++;
		occupationTotalHoles++;
	}
	
	public void remove_hole(Hole h, int band) {
		holesOnNP[band].remove(h);
		occupationVB[band]--;
		occupationTotalHoles--;
		
	}
	
	public void add_hop(Nanoparticle j) {
		Integer oldHop = electronHopMap.get(j);
		if( oldHop != null) {
			Integer newHop = oldHop + 1;
			electronHopMap.replace(j, newHop);
		}
		else {
			electronHopMap.put(j, 1);
		}
	}
	
	
    // update events, single NP version
    public void updateEvents(boolean self, Sample sample) {
    	// update self events
    	//System.out.println("updating NP"+ nanop);
    	if(self){
			for(int band=0; band<sample.nbnd; band++){
				for(Electron e: electronsOnNP[band]){
					// Erase all existing rates on the charge
					sample.rateOnSample -= e.ratesOnCharge;
					// look for new events and update the system total rate
					sample.rateOnSample += e.lookForEvents(sample);
				}
			}
		}
		// update neighboring events
		for(Nanoparticle neighborNP: nearestNeighbors){
			//// only update the ones that has not been updated
			////if(!currentUpdate.contains(neighborNP)){
				for(int band=0; band<sample.nbnd; band++){
					for(Electron e: neighborNP.electronsOnNP[band]){
						// Same as previous
						sample.rateOnSample -= e.ratesOnCharge;
						sample.rateOnSample += e.lookForEvents(sample);
					}
				}
				// add the updated NP to the updated list
				////currentUpdate.add(neighborNP);
			////}
		}
	}
    
	
	
	
	
}
