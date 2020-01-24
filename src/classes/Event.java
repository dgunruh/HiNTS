package classes;

public class Event {
	
    //int electronindex;
    //int holeindex  ;
    //int chargeindex ;    
    //int sourceparticle;          
    //int sourceorbital  ;        
    //int targetparticle;
    //int targetorbital  ;
    double energydiff;
    double rate;
    double poisson;
    double electronhole;
    int hopping_type;
    //double totalrate;
    int intra;
    double gap;
    
    public Event copy()  {
		try {
			return (Event) this.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
}
