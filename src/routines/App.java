package routines;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;


import classes.Sample;
import util.Configuration;
import util.Utility;

class Processor implements Callable<Double[]> {
	private int id;
	private int order_id;
	private boolean reverse_the_sample;
	private double temp;
	private int nelec;
	private double screeningFactor;
	private boolean randomSeed;
	private boolean necking;
	private double voltageRatio = 1.0;
	private String folderName;
	Double[] result;
	
	//here we allow different calls by overloading the constructor
	//base processor
	public Processor(String folderName, int id, boolean rSample, double temp, int nelec, double screeningFactor, boolean randomSeed, boolean necking){
		this.id = id;
		this.order_id = id;
		this.reverse_the_sample = rSample;
		this.temp = temp;
		this.nelec = nelec;
		this.screeningFactor = screeningFactor;
		this.randomSeed = randomSeed;
		this.necking = necking;
		this.folderName = folderName;
	}
	
	//voltage processor
	public Processor(String folderName, int id, boolean rSample, double temp, int nelec, double screeningFactor, boolean randomSeed, boolean necking, double voltageRatio){
		this.id = id;
		this.order_id = id;
		this.reverse_the_sample = rSample;
		this.temp = temp;
		this.nelec = nelec;
		this.screeningFactor = screeningFactor;
		this.randomSeed = randomSeed;
		this.necking = necking;
		this.voltageRatio = voltageRatio; //this is the unique feature
		this.folderName = folderName;
		
	}
	
	//repeat sample simulations processor
	public Processor(String folderName, int id, int order_id, boolean rSample, double temp, int nelec, double screeningFactor, boolean randomSeed, boolean necking){
		this.id = id;
		this.order_id = order_id; //this is the unique feature
		this.reverse_the_sample = rSample;
		this.temp = temp;
		this.nelec = nelec;
		this.screeningFactor = screeningFactor;
		this.randomSeed = randomSeed;
		this.necking = necking;
		this.folderName = folderName;
	}

	@Override
	public Double[] call() throws FileNotFoundException {
		System.out.println("Starting: "+id);
		
		result = new Double[2];
		
		Map<String, Object> params = new HashMap<>();
		
		int nholes = 0; //(int) Math.round(nelec/4.0); for troubleshooting purposes
		params.put("nelec", nelec);
		params.put("nholes", nholes);
		params.put("expected_nnanops", 800); //was 1200
		params.put("folderName", folderName);
		params.put("sample_no", id);
		params.put("reverse_the_sample", reverse_the_sample);
		params.put("feature", "mobility");
		params.put("voltage_ratio", voltageRatio);
		params.put("temperature", temp);
		params.put("closeNeighbor_thr", 0.0);
		params.put("neckedNeighbor_thr", 0.7); //0.42 allows the O_E to go up to kT. //Was 0.7
		params.put("screeningFactor", screeningFactor);
		params.put("randomSeed", randomSeed);
		params.put("necking", necking);
		
		Sample newsample = new Sample(params);
		
		// the first element is the id, the second element is the actual result
		result[0] = (double) order_id;
		result[1] = newsample.simulation();
		
		System.out.println("Completed: "+id);	
		
		return result;
	}
}

public class App {

	public static void main(String[] args) throws FileNotFoundException {
		
		int numberOfSamples = 200;
		String numberOfNanoparticles = "800"; //expected number
		int sampleStart = 0;
		String npArraySize = "20x2x20";
		String eDensity = ".00201";
		String diameter = "6.5";
		String defect = "twins";
		double T = 100.0;
		boolean randomSeed = false; //whether our MersenneTwisterFast will always use the same seeds, or random seeds
		boolean necking = false; //whether the simulation samples will have necking or not
		boolean preFlippedSamples = true;
		
		//These are the various types of runs that can be conducted
		boolean multipleDiameters = false; //cycle through batches of NPs at different diameters
		boolean repeatCalculations = false; //repeat mobility calculations for the same samples. Use random number as seed
		boolean linearElectronFilling = false; //cycle through electron filling for samples
		boolean ivCurve = false; //generate a curve of current versus voltage
		boolean multipleTemperatures = true; //cycle through multiple temperatures
		//if none of the above are "true", then will simply go through single batch of samples with no frills
		
		if(multipleDiameters) {
			List<Integer> nelecList = Arrays.asList(76, 109, 150, 199, 259, 329, 411, 506, 614, 736, 874, 1028, 1200, 1389, 1597, 1825, 2073); //.0015e-/nm^3 .5 ll
			
			for(double screeningFactor = 1.0; screeningFactor <= 1.0; screeningFactor += 1.0) {
				int j = 12;
				
				String title = "Cubic Results" + "\\" + Configuration.latticeStructure + "WITHHOLESnew_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_screeningFactor" + screeningFactor +".txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("Diameter(nm) Mobility(cm^2/Vs) Pair_STD");
				
				for(double diam = 9; diam <= 11; diam = diam + 0.5) {
					System.out.println("nanoparticle diameter: " + String.valueOf(diam));
					double nelec_doub = nelecList.get(j)*.0016/.0015;
					int nelec = (int) Math.round(nelec_doub);
					System.out.println("Number of electrons is " + nelec);
					ExecutorService executor = Executors.newFixedThreadPool(4);
					List<Future<Double[]>> futures = new ArrayList<>();
					
					String folderName = numberOfNanoparticles +"_"+diam + "nm";
					for(int i=0; i<numberOfSamples; i++){
						if(preFlippedSamples || (!preFlippedSamples && i%2==0)) {
							futures.add(executor.submit(new Processor(folderName, i,false, T,nelec,screeningFactor, randomSeed, necking)));
						}
						else {
							futures.add(executor.submit(new Processor(folderName, i,true, T,nelec,screeningFactor, randomSeed, necking)));
						}
						
					}
		
					// Stop accepting new tasks. Wait for all threads to terminate.
					executor.shutdown();
		
					System.out.println("All tasks submitted.");
		
					// wait for the processes to finish. Setting a time out limit of 1 day.
					try {
						executor.awaitTermination(1, TimeUnit.DAYS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
		
		
					System.out.println("All tasks completed.");
		
					List<Double[]> resultRaw = new ArrayList<>();
		
					for(Future<Double[]> future: futures){
						try {
							resultRaw.add(future.get());
						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}
					}
		
					// Average over regular and inverted samples
					System.out.println(resultRaw);
					double[] resultProcessed = Utility.processResult(resultRaw);
		
					double total = 0.0;
					for(int i = 0; i < resultProcessed.length; i++) {
						total += resultProcessed[i];
					}
					double mean = total/resultProcessed.length;
					
					double variance = 0.0;
					for(int i=0; i< resultProcessed.length;i++) {
						variance += Math.pow((resultProcessed[i]-mean),2);
					}
					double std = Math.sqrt(variance/(resultProcessed.length-1));
					double std_mean =  std/Math.sqrt(resultProcessed.length);
		
					String toAdd = String.valueOf(diam) + " " + mean + " " + std_mean;
					writer.println(toAdd);
					j += 1;
				}
				writer.close();
			}
		}
		
		else if(linearElectronFilling) {
			double screeningFactor = 1.0;
			String title = "Cubic Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_screeningFactor" + screeningFactor +".txt";
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs)");
			
			for(double fraction = 0.8; fraction <= 1.0; fraction = fraction + 0.1) {
				int nelec = (int) Math.round(fraction*800);
				int crack = 0;
				double diam = 6.5;
				System.out.println("Number Electrons: " + String.valueOf(nelec));
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
	
	
				String folderName = numberOfNanoparticles +"_"+diam + "nm";
				for(int i=0; i<numberOfSamples; i++){
					if(preFlippedSamples || (!preFlippedSamples && i%2==0)) {
						futures.add(executor.submit(new Processor(folderName, i,false,T,nelec,screeningFactor, randomSeed, necking)));
					}
					else {
						futures.add(executor.submit(new Processor(folderName, i,true,T,nelec,screeningFactor, randomSeed, necking)));
					}
				}
	
				// Stop accepting new tasks. Wait for all threads to terminate.
				executor.shutdown();
	
				System.out.println("All tasks submitted.");
	
				// wait for the processes to finish. Setting a time out limit of 1 day.
				try {
					executor.awaitTermination(1, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	
	
				System.out.println("All tasks completed.");
	
				List<Double[]> resultRaw = new ArrayList<>();
	
				for(Future<Double[]> future: futures){
					try {
						resultRaw.add(future.get());
					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}
	
				// Average over regular and inverted samples
				System.out.println(resultRaw);
				double[] resultProcessed = Utility.processResult(resultRaw);
				double total = 0.0;
				for(int i = 0; i < resultProcessed.length; i++) {
					total += resultProcessed[i];
					//String toAdd = String.valueOf(i) + " " + resultProcessed[i];
					//writer.println(toAdd);
				}
				double average = total/resultProcessed.length;
				String toAdd = String.valueOf(nelec) + " " + average;
				writer.println(toAdd);
			}
			writer.close();
		}
		
		else if(multipleTemperatures) {
			double screeningFactor = 1.0;
			List<String> disorderList = Arrays.asList(".00D",".015D", ".025D");
			for(String disorder: disorderList) {
				//String disorder = ".01D";
				String title = "TriDENS Results" + "/" + Configuration.latticeStructure + "_" + diameter + "nm" +npArraySize + "_" + eDensity + "_"  + "disorder" + disorder +".txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("Temperature(K) Mobility(cm^2/Vs) STD");
				
				for(double temp = 65.0; temp <= 80.0; temp += 2.5) {
					int nelec = 801; //1 e- per NP
					int crack = 0;
					double diam = 6.5;
					System.out.println("Temperature: " + temp);
					ExecutorService executor = Executors.newFixedThreadPool(4);
					List<Future<Double[]>> futures = new ArrayList<>();
		
		
					String folderName = numberOfNanoparticles +"_"+diam + "nm_200N_" + disorder;
					for(int i=0; i<numberOfSamples; i++){
						if(preFlippedSamples || (!preFlippedSamples && i%2==0)) {
							futures.add(executor.submit(new Processor(folderName, i,false,temp,nelec,screeningFactor, randomSeed, necking)));
						}
						else {
							futures.add(executor.submit(new Processor(folderName, i,true,temp,nelec,screeningFactor, randomSeed, necking)));
						}
					}
		
					// Stop accepting new tasks. Wait for all threads to terminate.
					executor.shutdown();
		
					System.out.println("All tasks submitted.");
		
					// wait for the processes to finish. Setting a time out limit of 1 day.
					try {
						executor.awaitTermination(1, TimeUnit.DAYS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
		
		
					System.out.println("All tasks completed.");
		
					List<Double[]> resultRaw = new ArrayList<>();
		
					for(Future<Double[]> future: futures){
						try {
							resultRaw.add(future.get());
						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}
					}
		
					// Average over regular and inverted samples
					System.out.println(resultRaw);
					double[] resultProcessed = Utility.processResult(resultRaw);
					double total = 0.0;
					int nonZeroResults = 0;
					for(int i = 0; i < resultProcessed.length; i++) {
						if(resultProcessed[i] > 0) {
							total += resultProcessed[i];
							nonZeroResults++;
						}
						//total += resultProcessed[i];
						//String toAdd = String.valueOf(i) + " " + resultProcessed[i];
						//writer.println(toAdd);
					}
					//double average = total/resultProcessed.length;
					double average = total/nonZeroResults;
					
					double variance = 0.0;
					for(int i=0; i< nonZeroResults;i++) {
						variance += Math.pow((resultProcessed[i]-average),2);
					}
					double std = Math.sqrt(variance/(nonZeroResults-1));
					double std_mean =  std/Math.sqrt(nonZeroResults);
					String toAdd = String.valueOf(temp) + " " + average +" " + std_mean;
					writer.println(toAdd);
				}
				writer.close();
			}
		}
		
		else if(ivCurve) {
			double screeningFactor = 1.0;
			
			for(double voltageRatio = 0.25; voltageRatio <= 1.5; voltageRatio = voltageRatio + 0.25) {
				String title = "TriDENS Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_voltageRatio" + voltageRatio +".txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("SampleNumber Mobility(cm^2/Vs)");
				int nelec = 350;
				int crack = 0;
				double diam = 6.5;
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
	
	
				String folderName = numberOfNanoparticles +"_"+diam + "nm";
				for(int i=0; i<numberOfSamples; i++){
					if(preFlippedSamples || (!preFlippedSamples && i%2==0)) {
						futures.add(executor.submit(new Processor(folderName, i,false,T,nelec,screeningFactor, randomSeed, necking, voltageRatio)));
					}
					else {
						futures.add(executor.submit(new Processor(folderName, i,true,T,nelec,screeningFactor, randomSeed, necking, voltageRatio)));
					}
				}
	
				// Stop accepting new tasks. Wait for all threads to terminate.
				executor.shutdown();
	
				System.out.println("All tasks submitted.");
	
				// wait for the processes to finish. Setting a time out limit of 1 day.
				try {
					executor.awaitTermination(1, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	
	
				System.out.println("All tasks completed.");
	
				List<Double[]> resultRaw = new ArrayList<>();
	
				for(Future<Double[]> future: futures){
					try {
						resultRaw.add(future.get());
					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}
	
				// Average over regular and inverted samples
				System.out.println(resultRaw);
				double[] resultProcessed = Utility.processResult(resultRaw);

				for(int i = 0; i < resultProcessed.length; i++) {
					String toAdd = String.valueOf(i) + " " + resultProcessed[i];
					writer.println(toAdd);
				}
				writer.close();
			}
		}
		
		//in this scenario, we will assume that we are examining individual samples, so repeat calculations several times for each
		else if(repeatCalculations) {
			int nelec = 200; //.002 e-/nm^3 triclinic
			double diam = 6.5;
			
			int repetitions = 10;
			double screeningFactor = 1.0;
			//String title = "Cubic Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_screeningFactor" + screeningFactor +".txt";
			//String title = "Cubic Results" + "\\" + "MouleMobilityResult.txt";
			String title = "TriDENS Results" + "\\" + Configuration.latticeStructure + "_" + diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K" + defect +".txt";
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs) Pair_STD");
			
			for(int k = 0; k < 1; k++) {
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
				String folderName = numberOfNanoparticles +"_"+diam + "nm";
				
				//submit jobs with an ordering id separate from their actual sample id
				//this will allow us to average over groups of runs on the same sample
				for(int i=0; i<numberOfSamples; i+=2){
					for(int j=0; j<repetitions; j++) {
						int nonflippedID = i*repetitions + 2*j;
						int flippedID = i*repetitions + 2*j + 1;
						
						if(preFlippedSamples) {
							futures.add(executor.submit(new Processor(folderName, i,nonflippedID, false, T,nelec,screeningFactor, randomSeed, necking)));
							futures.add(executor.submit(new Processor(folderName, i+1,flippedID, false, T,nelec, screeningFactor, randomSeed, necking)));
						}
						else {
							futures.add(executor.submit(new Processor(folderName, i,nonflippedID, false, T,nelec,screeningFactor, randomSeed, necking)));
							futures.add(executor.submit(new Processor(folderName, i+1,flippedID, true, T,nelec, screeningFactor, randomSeed, necking)));
						}
					}
				}
	
				// Stop accepting new tasks. Wait for all threads to terminate.
				executor.shutdown();
	
				System.out.println("All tasks submitted.");
	
				// wait for the processes to finish. Setting a time out limit of 1 day.
				try {
					executor.awaitTermination(1, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	
	
				System.out.println("All tasks completed.");
	
				List<Double[]> resultRaw = new ArrayList<>();
	
				for(Future<Double[]> future: futures){
					try {
						resultRaw.add(future.get());
					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}
				
				for(int i = 0; i < resultRaw.size(); i++) {
					System.out.println("Result: " + resultRaw.get(i));
				}
	
				// Average over regular and inverted samples
				double[] resultProcessed = Utility.processResult(resultRaw);
				System.out.println("Result processed length: " + resultProcessed.length);
	
				double[] total = new double[resultProcessed.length/repetitions];
				double[] mean = new double[total.length];
				
				for(int i = 0; i < resultProcessed.length; i++) {
					System.out.println("Result: " + resultProcessed[i]);
				}
				
				for(int i = 0; i < resultProcessed.length/repetitions; i++) {
					for(int j=0; j < repetitions; j++) {
						total[i] += resultProcessed[i*repetitions + j];
						System.out.println("Total is: " + total[i]);
					}
					mean[i] = total[i]/repetitions;
					System.out.println("Mean is: " + mean[i]);
				}
				
				double[] variance = new double[mean.length];
				double[] std_mean = new double[mean.length];
				
				for(int i=0; i< mean.length;i++) {
					for(int j = 0; j<repetitions; j++) {
						variance[i] += Math.pow((resultProcessed[i*repetitions + j]-mean[i]),2);
						System.out.println("Variance is: " + variance[i]);
					}
					std_mean[i] = Math.sqrt(variance[i]/(repetitions-1));
					std_mean[i] = std_mean[i]/Math.sqrt(repetitions);
				}
				
				for(int i = 0; i < mean.length; i++) {
					double sampleResult = mean[i];
					double sampleSTD = std_mean[i];
					String toAdd = i + " " + sampleResult + " " + sampleSTD;
					writer.println(toAdd);
				}
			}
			writer.close();
		}
		
		else {
			int nelec = 300; //.002 e-/nm^3 triclinic 20x2x20
			int nhole = 0;
			double diam = 6.0;
			double screeningFactor = 1.0;
			
			String title = "TriDENS Results" + "\\" + Configuration.latticeStructure + "_" + diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K" + numberOfSamples+ defect + "TEST"+ ".txt";
			
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs)");
			
			
			ExecutorService executor = Executors.newFixedThreadPool(4);
			List<Future<Double[]>> futures = new ArrayList<>();

			String folderName = numberOfNanoparticles +"_"+diam + "nm";
			for(int i=0; i<numberOfSamples; i++){
				if(preFlippedSamples || (!preFlippedSamples && i%2==0)) {
					futures.add(executor.submit(new Processor(folderName, i + sampleStart,i, false, T,nelec, screeningFactor, randomSeed, necking)));
				}
				else {
					futures.add(executor.submit(new Processor(folderName, i + sampleStart,i, true, T,nelec, screeningFactor, randomSeed, necking)));
				}
			}

			// Stop accepting new tasks. Wait for all threads to terminate.
			executor.shutdown();

			System.out.println("All tasks submitted.");

			// wait for the processes to finish. Setting a time out limit of 1 day.
			try {
				executor.awaitTermination(1, TimeUnit.DAYS);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}


			System.out.println("All tasks completed.");

			List<Double[]> resultRaw = new ArrayList<>();

			for(Future<Double[]> future: futures){
				try {
					resultRaw.add(future.get());
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}

			// Average over regular and inverted samples
			double[] resultProcessed = Utility.processResult(resultRaw);

			for(int i = 0; i < resultProcessed.length; i++) {
				double sampleResult = resultProcessed[i];
				String toAdd = i + " " + sampleResult;
				writer.println(toAdd);
			}
			writer.close();
		}
		
	}

}



