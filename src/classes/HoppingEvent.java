package classes;

public class HoppingEvent extends Event {
	
	Electron hostElectron;
	Hole hostHole;
	Nanoparticle targetNP, sourceNP;
	int targetOrbital, sourceOrbital;
	String type;
	boolean intra; // if intracluster hopping
	 
	public HoppingEvent () {
		type = "empty";
	}
	
	
	public HoppingEvent(Electron e, Nanoparticle target, Nanoparticle source, int to, int so, boolean intra){
		hostElectron = e;
		targetNP = target;
		sourceNP = source;
		targetOrbital = to;
		sourceOrbital = so;
		if(!intra)
			type = "regular electron hopping";
		else
			type = "cluster electron hopping";
	}
	
	public HoppingEvent(Hole h, Nanoparticle target, Nanoparticle source, int to, int so, boolean intra){
		hostHole = h;
		targetNP = target;
		sourceNP = source;
		targetOrbital = to;
		sourceOrbital = so;
		if(!intra)
			type = "regular hole hopping";
		else
			type = "cluster hole hopping";

	}
	
	public void setRate(double rate) {
		if(rate<0)
			System.out.println("rate cannot be negative");
		this.rate = rate;
	}
}
