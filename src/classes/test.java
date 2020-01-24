package classes;


public class test {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int i;
		
		i = 0;
		
		while(i<10){
			
			for(int j=0; j<50; j++){
				i += j;
				
				if(i>=10)
					System.out.println("already"+i+" "+j);
				if(j==49)
					System.out.println("full circle");
			}
		}
		
		System.out.print(i);

	}

}
