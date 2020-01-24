package classes;

import java.util.ArrayList;

import util.Configuration;
import util.Constants;

public class Hole extends Charge {
	public Hole() {
		charge = 1*Constants.sqrt2;
		mass = Configuration.hmass;
		ratesOnCharge = 0;
	}



	public double lookForEvents(Sample sample) {
		//System.out.println(this.hostNP);
		double degeneracy;
		double netDistance;
		double newRate;
		HoppingEvent newEvent;
		Nanoparticle sourceNP = this.hostNP;
		int sourceOrbital = this.hostOrbital;
		
		this.hoppings = new ArrayList<HoppingEvent>();
		
		// zero total rates on charge
		ratesOnCharge = 0 ;
		
		// add cluster hoppings first, if energy levels are close enough, then metallic/ ohmic transport
		for(Nanoparticle neighborNP : sourceNP.holeCloseNeighbors){
			for(int band=0; band<sample.nbnd; band++){
				// band is targetOrbital
				// allow hopping towards SOURCE only
				netDistance = (neighborNP.z-sourceNP.z-sample.cellz*Math.round((neighborNP.z-sourceNP.z)/sample.cellz));
				if(netDistance<0 && neighborNP.occupationVB[band]<neighborNP.occupationMAX[band] ){
					
					newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, true);
					degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationVB[band]) ; 
					newRate = charge*degeneracy*(sample.metalPrefactor)*sample.voltage/sample.cellz*netDistance;
					newEvent.setRate(newRate);
					this.hoppings.add(newEvent);
					ratesOnCharge += newRate ;
				}
			}
		}
		
		//add necking transport next
		if(sourceNP.neckedNeighbors.length != 0) {
			for(Nanoparticle neighborNP : sourceNP.neckedNeighbors) {
				for(int band = 0; band < sample.nbnd; band++) {
					//band is target orbital
					double energy_diff = 0.0;

					// Kinetic energy difference
					energy_diff += neighborNP.vbenergy[band] - sourceNP.vbenergy[sourceOrbital] ;

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

					if(energy_diff <= overlap_energy && neighborNP.occupationVB[band]<neighborNP.occupationMAX[band] ){

						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationVB[band]) ; 

						double thisNP_electron_density = sourceNP.occupationTotalElectron/(Math.PI/6.0*Math.pow(sourceNP.diameter, 3.0));
						double neighborNP_electron_density = (neighborNP.occupationTotalElectron + 1)/(Math.PI/6.0*Math.pow(neighborNP.diameter, 3.0));
						double electron_density = (thisNP_electron_density + neighborNP_electron_density)/2.0;
						double matrixElement = 9*electron_density*Math.pow(neckRadius, 3.0)/(Configuration.emass*Math.pow(diameter, 2.0));  //9*hbar^2*n*rho^3/(m^* * d^2)

						newRate = 2*Math.PI*Math.pow(matrixElement, 2.0)*degeneracy;
						newEvent.setRate(newRate);
						this.hoppings.add(newEvent);
						ratesOnCharge += newRate ;
					}
					else if(energy_diff > overlap_energy && neighborNP.occupationVB[band]<neighborNP.occupationMAX[band]) {
						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationVB[band]) ; 
						double thisNP_electron_density = sourceNP.occupationTotalElectron/(Math.PI/6.0*Math.pow(sourceNP.diameter, 3.0));
						double neighborNP_electron_density = (neighborNP.occupationTotalElectron + 1)/(Math.PI/6.0*Math.pow(neighborNP.diameter, 3.0));
						double electron_density = (thisNP_electron_density + neighborNP_electron_density)/2.0;
						double matrixElement = 9*electron_density*Math.pow(neckRadius, 3.0)/(Configuration.emass*Math.pow(diameter, 2.0));

						double thermal_activation = Math.exp(-energy_diff/sample.temperature);
						//System.out.println("Thermal activation is: " + thermal_activation);
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
			for(int band=0; band<sample.nbnd; band++){
				if(neighborNP.occupationVB[band]<neighborNP.occupationMAX[band]){
					newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
					newRate = calculateRate(neighborNP, sourceNP, band, sourceOrbital, sample);
					newEvent.setRate(newRate);
					this.hoppings.add(newEvent);
					ratesOnCharge += newRate;
				}
			}
		}
		return ratesOnCharge;
	}

	
	private double calculateCharging(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample){
	    /* OLD METHOD
		double charging = 0.0;
    	// On-Site charging energy
    	if(sample.lcapacitance0)
    		charging += targetNP.hselfenergy0 - sourceNP.hselfenergy0 ;
    	
    	charging += targetNP.occupationTotalHoles*targetNP.hselfenergy - (sourceNP.occupationTotalHoles-1)*sourceNP.hselfenergy;
    	return charging;
    	*/
		
		double chargingBefore = (targetNP.calculateOnSiteHoleCharging(targetNP.occupationTotalElectron)
                + sourceNP.calculateOnSiteHoleCharging((sourceNP.occupationTotalElectron)));
	    double chargingAfter = (targetNP.calculateOnSiteHoleCharging(targetNP.occupationTotalElectron+1) 
	    		+ sourceNP.calculateOnSiteHoleCharging((sourceNP.occupationTotalElectron-1)));
	    
	    return (chargingAfter-chargingBefore)/sample.screeningFactor;
    }
	
	private double calculateExciton(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample){
	    /* OLD METHOD
		double exciton = 0.0;
    	
	    double excitonSP = Math.min(sourceNP.occupationTotalHoles, sourceNP.occupationTotalElectron);
	    double excitonTP = Math.min(targetNP.occupationTotalHoles, targetNP.occupationTotalElectron+1);
	    
	    exciton = (excitonTP>0) ? -targetNP.hselfenergy0 : 0 ;
	    exciton += (excitonSP>0) ? sourceNP.hselfenergy0 : 0 ; 
	    
	    excitonSP = Math.min(sourceNP.occupationTotalHoles, sourceNP.occupationTotalElectron-1);
	    excitonTP = Math.min(targetNP.occupationTotalHoles, targetNP.occupationTotalElectron);

	    exciton += -excitonTP*targetNP.hselfenergy + excitonSP*sourceNP.hselfenergy ;
	    
    	return exciton;
    	*/
		
		//new method
	    double chargingBefore = (targetNP.calculateOnSiteElectronCharging(targetNP.occupationTotalElectron)
                + sourceNP.calculateOnSiteElectronCharging((sourceNP.occupationTotalElectron)));
	    double chargingAfter = (targetNP.calculateOnSiteElectronCharging(targetNP.occupationTotalElectron+1) 
	    		+ sourceNP.calculateOnSiteElectronCharging((sourceNP.occupationTotalElectron-1)));
	    
	    return (chargingAfter-chargingBefore)/sample.screeningFactor;
    }

	
	// !! UNRELIABLE PART !!
	private double nearestNeighborPoisson(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample) {
		double nnPoisson=0;
		double ccdistance; 
	    // for the electron before hopping this comes with negative sign as this is subtracted
		for(Nanoparticle sourceNeighbor: sourceNP.nearestNeighbors){
			ccdistance = sourceNP.centerDistanceMap.get(sourceNeighbor);
			nnPoisson += -Constants.e2*(sourceNeighbor.occupationTotalElectron-sourceNeighbor.occupationTotalHoles) / sample.dcout / ccdistance;
		}
		// for the electron after hopping
		// CAUTION! here should be for the hole after the hopping !!
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
    	// Kinetic energy difference
    	energy_diff += targetNP.vbenergy[targetOrbital] - sourceNP.vbenergy[sourceOrbital] ;
    	
    	// On-Site charging energy
    	energy_diff += calculateCharging(targetNP, sourceNP, sample);
    	
    	// Electron-hole exciton interaction
    	energy_diff += calculateExciton(targetNP, sourceNP, sample);
    	
    	// External potential
        energy_diff += charge*(sample.voltage/sample.cellz)*(targetNP.z-sourceNP.z-sample.cellz*Math.round((targetNP.z-sourceNP.z)/sample.cellz));

    	// Nearest neighbor interaction
        if(sample.poissonSolver=="nn")
        	energy_diff += nearestNeighborPoisson(targetNP, sourceNP, sample);
        
        //question: should this be vb energy, or should it be cb energy? Luman had it as cb energy
        overlap = Math.sqrt(-sample.emass*(sourceNP.vbenergy[sourceOrbital] + targetNP.vbenergy[targetOrbital]) / 2.0) ;
        
        energy_diff = Math.abs(energy_diff); //since energy is negative, make it the absolute value

        // Two cases: Miller-Abrahms or Marcus
        switch (sample.hoppingMechanism) {
		case "ma":
			if (npdistance > 0) { rate = sample.jumpfreq * Math.exp(-2.0*npdistance*overlap); 
			}
			else { 
				rate = sample.jumpfreq;
				//rate = 0.0;
				//System.out.println("hitting it");
			}
			
			if(energy_diff>0.0)
				rate = rate*Math.exp(-energy_diff/sample.temperature);
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
    	sourceNP.remove_hole(this, sourceBand);
    	// add to the new NP
    	this.hostNP = destinationNP;
    	this.hostOrbital = destinationBand;
    	destinationNP.add_hole(this, destinationBand);

    	double[] sourceDestination = {destinationNP.x, destinationNP.z};
    	this.visitedNPs.add(sourceDestination);
    }

}

