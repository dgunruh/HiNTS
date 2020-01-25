package classes;

import java.util.ArrayList;
import java.util.Arrays;

import util.Configuration;
import util.Constants;

public class Electron extends Charge {
	public Electron(){
		charge = -1.0*Constants.sqrt2;
		mass = Configuration.emass;
		ratesOnCharge = 0.0;
		hoppings = new ArrayList<HoppingEvent>();
		visitedNPs = new ArrayList<double[]>();
		//events = new RegularHopping[5];
	}
	
	
	
	public double lookForEvents(Sample sample) {
		double degeneracy;
		double netDistance;
		double newRate;
		HoppingEvent newEvent;
		Nanoparticle sourceNP = this.hostNP;
		int sourceOrbital = this.hostOrbital;
		
		this.hoppings.clear();
		
		// zero total rates on charge
		ratesOnCharge = 0 ;

		// add cluster hoppings next
		for(Nanoparticle neighborNP : sourceNP.closeNeighbors){
			for(int band=0; band<sample.nbnd; band++){
				// band is targetOrbital
				// allow hopping towards drain only
				netDistance = (neighborNP.z-sourceNP.z-sample.cellz*Math.round((neighborNP.z-sourceNP.z)/sample.cellz));
				if(netDistance>0 && neighborNP.occupationCB[band]<neighborNP.occupationMAX[band] ){
					newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, true);
					degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]) ; 
					newRate = -charge*degeneracy*(sample.metalPrefactor)*sample.voltage/sample.cellz*netDistance;
					newEvent.setRate(newRate);
					this.hoppings.add(newEvent);
					ratesOnCharge += newRate ;
					System.out.println("CLUSTER RATE IS: " + newRate);
				}
			}
		}
		
		//add necking transport next
		if(sourceNP.neckedNeighbors.length != 0) {
			System.out.println("Adding necking");
			for(Nanoparticle neighborNP : sourceNP.neckedNeighbors) {
				System.out.println("Necked transport!");
				for(int band = 0; band < sample.nbnd; band++) {
					//band is target orbital
					double energy_diff = 0.0;
					
					// Kinetic energy difference
			    	energy_diff += neighborNP.cbenergy[band] - sourceNP.cbenergy[sourceOrbital] ;
			    	
			    	// On-Site charging energy
			    	energy_diff += calculateCharging(neighborNP, sourceNP, sample);
			    	
			    	// Electron-hole exciton interaction
			    	energy_diff += calculateExciton(neighborNP, sourceNP, sample);
	
					// allow hopping only towards nanoparticles with a lower energy level
					netDistance = (neighborNP.z-sourceNP.z-sample.cellz*Math.round((neighborNP.z-sourceNP.z)/sample.cellz));
					energy_diff += charge*sample.voltage*netDistance/sample.cellz;
					
					double neckRadius = sourceNP.neckRadiusMap.get(neighborNP);
					double diameter = (sourceNP.diameter + neighborNP.diameter)/2.0;
					double overlap_energy = sample.ediff_neck_thr*2*neckRadius/diameter; //Scale proportionally to neck thickness
					//System.out.println("Overlap energy is: " + overlap_energy);
					//System.out.println("kT is: " + sample.temperature);
					//double overlap_energy = sample.temperature
					
					if(energy_diff <= overlap_energy && neighborNP.occupationCB[band]<neighborNP.occupationMAX[band] ){
						
						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]) ; 
						
						double thisNP_electron_density = sourceNP.occupationTotalElectron/(Math.PI/6.0*Math.pow(sourceNP.diameter, 3.0));
						double neighborNP_electron_density = (neighborNP.occupationTotalElectron + 1)/(Math.PI/6.0*Math.pow(neighborNP.diameter, 3.0));
						//double electron_density = sample.nelec/(sample.cellx*sample.celly*sample.cellz);
						double electron_density = (thisNP_electron_density + neighborNP_electron_density)/2.0;
						double matrixElement = 9*electron_density*Math.pow(neckRadius, 3.0)/(Configuration.emass*Math.pow(diameter, 2.0));  //9*hbar^2*n*rho^3/(m^* * d^2)
						
						newRate = 2*Math.PI*Math.pow(matrixElement, 2.0)*degeneracy;
						newEvent.setRate(newRate);
						this.hoppings.add(newEvent);
						ratesOnCharge += newRate ;
					}
					else if(energy_diff > overlap_energy && neighborNP.occupationCB[band]<neighborNP.occupationMAX[band]) {
						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]) ; 
						double thisNP_electron_density = sourceNP.occupationTotalElectron/(Math.PI/6.0*Math.pow(sourceNP.diameter, 3.0));
						double neighborNP_electron_density = (neighborNP.occupationTotalElectron + 1)/(Math.PI/6.0*Math.pow(neighborNP.diameter, 3.0));
						//double electron_density = sample.nelec/(sample.cellx*sample.celly*sample.cellz);
						double electron_density = (thisNP_electron_density + neighborNP_electron_density)/2.0;
						double matrixElement = 9*electron_density*Math.pow(neckRadius, 3.0)/(Configuration.emass*Math.pow(diameter, 2.0));
						
						double thermal_activation = Math.exp(-energy_diff/sample.temperature);
						newRate = 2*Math.PI*Math.pow(matrixElement, 2.0)*degeneracy*thermal_activation; //make sure rate crosses over from hopping
						newEvent.setRate(newRate);
						this.hoppings.add(newEvent);
						ratesOnCharge += newRate ;
					}
				}
			}
		}
			
		// add regular hoppings last
		for(Nanoparticle neighborNP : sourceNP.nearestNeighbors){
			if(!Arrays.asList(sourceNP.neckedNeighbors).contains(neighborNP)) {
				for(int band=0; band<sample.nbnd; band++){
					if(neighborNP.occupationCB[band]<neighborNP.occupationMAX[band]){
						//HoppingEvent(Electron e, Nanoparticle target, Nanoparticle source, int to, int so, boolean intra)
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]);
						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						newRate = degeneracy*calculateRate(neighborNP, sourceNP, band, sourceOrbital, sample);
						newEvent.setRate(newRate);
						this.hoppings.add(newEvent);
						ratesOnCharge += newRate;
					}
				}
			}
			else {
				System.out.println("Contains neighbor!");
			}
		}
		
		return ratesOnCharge;
	}

	
	private double calculateCharging(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample){
    	// On-Site charging energy
		/*OLD METHOD
		 * 
		 * double charging = 0.0;
	    
    	if(sample.lcapacitance0)
    		charging += targetNP.selfenergy0 - sourceNP.selfenergy0 ;
    	
    	charging += targetNP.occupationTotalElectron*targetNP.selfenergy - (sourceNP.occupationTotalElectron-1)*sourceNP.selfenergy;
    	
	    */
		
		//new method
	    double chargingBefore = (targetNP.calculateOnSiteElectronCharging(targetNP.occupationTotalElectron)
                + sourceNP.calculateOnSiteElectronCharging((sourceNP.occupationTotalElectron)));
	    double chargingAfter = (targetNP.calculateOnSiteElectronCharging(targetNP.occupationTotalElectron+1) 
	    		+ sourceNP.calculateOnSiteElectronCharging((sourceNP.occupationTotalElectron-1)));
	    
	    return (chargingAfter-chargingBefore)/sample.screeningFactor;
    }
	
	private double calculateExciton(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample){
	    /*old method for calculating excitonic interactions
	     * 
	     * double exciton = 0.0;
    	
	    double excitonSP = Math.min(sourceNP.occupationTotalHoles, sourceNP.occupationTotalElectron);
	    double excitonTP = Math.min(targetNP.occupationTotalHoles, targetNP.occupationTotalElectron+1);
	    
	    exciton = (excitonTP>0) ? -targetNP.selfenergy0 : 0 ;
	    exciton += (excitonSP>0) ? sourceNP.selfenergy0 : 0 ; 
	    
	    excitonSP = Math.min(sourceNP.occupationTotalHoles, sourceNP.occupationTotalElectron-1);
	    excitonTP = Math.min(targetNP.occupationTotalHoles, targetNP.occupationTotalElectron);

	    exciton += -excitonTP*targetNP.selfenergy + excitonSP*sourceNP.selfenergy ;
	    
    	return exciton;*/
		
		//new method for calculating excitonic interactions
		//this uses the binding energy for excitons constrained to NPs, calculated under the effective mass approximation
		double excitonBefore = -targetNP.calculateOnSiteExcitonCharging(targetNP.occupationTotalElectron, targetNP.occupationTotalHoles, sample)
				- sourceNP.calculateOnSiteExcitonCharging(sourceNP.occupationTotalElectron, sourceNP.occupationTotalHoles, sample);
		double excitonAfter = -targetNP.calculateOnSiteExcitonCharging(targetNP.occupationTotalElectron + 1, targetNP.occupationTotalHoles, sample)
				- sourceNP.calculateOnSiteExcitonCharging(sourceNP.occupationTotalElectron - 1, sourceNP.occupationTotalHoles, sample);
		
		return (excitonAfter - excitonBefore)/sample.screeningFactor;
    }

	private double nearestNeighborPoisson(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample) {
		double nnPoisson=0;
		double ccdistance; 
	    // for the electron before hopping this comes with negative sign as this is subtracted
		for(Nanoparticle sourceNeighbor: sourceNP.nearestNeighbors){
			ccdistance = sourceNP.centerDistanceMap.get(sourceNeighbor);
			nnPoisson += -Constants.e2*(sourceNeighbor.occupationTotalElectron-sourceNeighbor.occupationTotalHoles) / sample.dcout / ccdistance;
		}
		// for the electron after hopping
		for(Nanoparticle targetNeighbor: targetNP.nearestNeighbors){
			ccdistance = targetNP.centerDistanceMap.get(targetNeighbor);
			
			if(targetNeighbor==sourceNP)
				nnPoisson += Constants.e2*(targetNeighbor.occupationTotalElectron - targetNeighbor.occupationTotalHoles -1) / sample.dcout / ccdistance;	
			else
				nnPoisson += Constants.e2*(targetNeighbor.occupationTotalElectron- targetNeighbor.occupationTotalHoles) / sample.dcout / ccdistance;
		}
		return nnPoisson;
	}

    private double calculateRate(Nanoparticle targetNP, Nanoparticle sourceNP, int targetOrbital, int sourceOrbital, Sample sample){
    	double energy_diff = 0.0;
    	double rate = 0, overlap;
    	double npdistance = targetNP.edgeDistanceMap.get(sourceNP);
    	if (npdistance < 0) {
    		System.out.println("nanoparticles too close");
    		System.out.println("NP 1: " + targetNP.id +" NP 2: " + sourceNP.id);
    		System.out.println("deltax: " + (targetNP.x - sourceNP.x)*Constants.bohrtonm);
    		System.out.println("deltay: " + (targetNP.y - sourceNP.y)*Constants.bohrtonm);
    		System.out.println("deltaz: " + (targetNP.z - sourceNP.z)*Constants.bohrtonm);
    		System.out.println(Math.min(Math.min(Math.abs(targetNP.z-sourceNP.z), Math.abs(targetNP.z-(sourceNP.z-sample.cellz))), Math.abs(targetNP.z-(sourceNP.z+sample.cellz)))*Constants.bohrtonm);
    		System.out.println("diameter 1: " + targetNP.diameter*Constants.bohrtonm + " diameter 2: " + sourceNP.diameter*Constants.bohrtonm);
    	}
    	//if (npdistance < 0) npdistance = Configuration.ligandLength*Constants.nmtobohr; 
    	// Kinetic energy difference
    	energy_diff += targetNP.cbenergy[targetOrbital] - sourceNP.cbenergy[sourceOrbital] ;
    	
    	// On-Site charging energy
    	energy_diff += calculateCharging(targetNP, sourceNP, sample);
    	//System.out.println("Charging energy is: " + calculateCharging(targetNP, sourceNP, sample)*Constants.rytoev);
    	
    	// Electron-hole exciton interaction
    	energy_diff += calculateExciton(targetNP, sourceNP, sample);
    	if(calculateExciton(targetNP, sourceNP, sample) != 0.0) System.out.println("exciton was used");
    	
    	// External potential
        energy_diff += charge*(sample.voltage/sample.cellz)*(targetNP.z-sourceNP.z-sample.cellz*Math.round((targetNP.z-sourceNP.z)/sample.cellz));

        //System.out.println("Energy diff is: " + energy_diff*Constants.rytoev);
        
    	// Nearest neighbor interaction
        if(sample.poissonSolver=="nn")
        	energy_diff += nearestNeighborPoisson(targetNP, sourceNP, sample);
        
        overlap = Math.sqrt(-sample.emass*(sourceNP.cbenergy[sourceOrbital] + targetNP.cbenergy[targetOrbital]) / 2.0) ;
        
        //System.out.println("Energy diff is: " + energy_diff*Constants.rytoev);
        //npdistance = .4*Constants.nmtobohr;
        //System.out.println("Overlap is: " + overlap*npdistance);
        //System.out.println("Source NP diameter:" + sourceNP.diameter*Constants.bohrtonm);

        // Two cases: Miller-Abrahms or Marcus
        switch (sample.hoppingMechanism) {
		case "ma":
			if (npdistance > 0) { rate = sample.jumpfreq * Math.exp(-2.0*npdistance*overlap); 
			}
			else { 
				rate = sample.jumpfreq;
				//rate = 0.0;
				System.out.println("hitting it");
			}
			
			//This is an implementation of the overlap energy for samples with disorder. A stopgap measure
			if(energy_diff > sample.ediff_hopping_thr) {
				rate = rate*Math.exp(-energy_diff/sample.temperature);
				//System.out.println("Rate is: " + rate);
				//System.out.println("temperature limited step!");
			}
			else {
				/*
				System.out.println("kT: " + sample.temperature*Constants.rytoev);
				System.out.println("NP separation: " + (targetNP.z - sourceNP.z)*Constants.bohrtonm);
				System.out.println("cb energy diff: " + (targetNP.cbenergy[targetOrbital] - sourceNP.cbenergy[sourceOrbital])*Constants.rytoev);
				System.out.println("charging energy diff: " + calculateCharging(targetNP, sourceNP, sample)*Constants.rytoev);
				System.out.println("voltage diff: " + charge*(sample.voltage/sample.cellz)*(targetNP.z-sourceNP.z-sample.cellz*Math.round((targetNP.z-sourceNP.z)/sample.cellz))*Constants.rytoev);
				System.out.println("WKB term: " + Math.exp(-2.0*npdistance*overlap));
				*/
				//only allow forward or sideways transport
				double netDistance = (targetNP.z-sourceNP.z-sample.cellz*Math.round((targetNP.z-sourceNP.z)/sample.cellz));
				if(netDistance < 0) {
					rate = 0.0;
				}
			}
			break;

		case "marcus":
			rate = Constants.tpi*sample.marcusprefac2*Math.exp(-2.0 * npdistance * overlap) / Math.sqrt(Constants.fpi*sample.reorgenergy*sample.temperature);
	        rate = rate*Math.exp(-Math.pow((sample.reorgenergy + energy_diff), 2.0) / (4.0*sample.reorgenergy*sample.temperature));
			break;
		}
    	
    	return rate;
    }
    
    
    public void move(Nanoparticle sourceNP, int sourceBand, Nanoparticle destinationNP, int destinationBand) {
		// remove from old NP first
    	sourceNP.remove_electron(this, sourceBand);
    	// add to the new NP
    	this.hostNP = destinationNP;
    	this.hostOrbital = destinationBand;
    	destinationNP.add_electron(this, destinationBand);
    	
    	double[] sourceDestination = {destinationNP.x, destinationNP.z};
    	this.visitedNPs.add(sourceDestination);
	}
    
    
    public double calculateEnergy() {
    	return 0;
		
	}


}
