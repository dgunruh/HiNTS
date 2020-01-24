package classes;

import java.util.ArrayList;
import java.util.List;

public class Charge {
    
	Nanoparticle hostNP;
    int hostOrbital;
    
    double charge;
    double mass;
    double ratesOnCharge;
    Event[] events;
    public ArrayList<HoppingEvent> hoppings;
    public ArrayList<double[]> visitedNPs;
    
    
    public void setHost(Nanoparticle host, int orbital){
    	hostNP = host;
    	hostOrbital = orbital;
    }
    
    public Charge clone() {
    	return this;
		
	}

}

