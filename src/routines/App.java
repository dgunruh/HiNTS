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
	public Processor(String folderName, int id, int order_id, boolean rSample, double temp, int nelec, double screeningFactor, boolean randomSeed, boolean necking, double voltageRatio){
		this.id = id;
		this.order_id = order_id;
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
		
		int numberOfSamples = 200; //total number of samples including flipped ones. If not pre-flipped, this is the total expected number
		String numberOfNanoparticles = "450"; //expected number
		int sampleStart = 0;
		String npArraySize = "20x2x20";
		String eDensity = ".0016";
		String diameter = "3.5";
		String defect = "twins";
		double T = 300.0;
		boolean randomSeed = false; //whether our MersenneTwisterFast will always use the same seeds, or random seeds
		boolean necking = false; //whether the simulation samples will have necking or not
		boolean preFlippedSamples = true; //whether pre-flipped pairs of samples are present or not. If not, do the flipping internally. 
		
		//These are the various types of runs that can be conducted
		boolean multipleDiameters = false; //cycle through batches of NPs at different diameters
		boolean repeatCalculations = false; //repeat mobility calculations for the same samples. Use random number as seed
		boolean linearElectronFilling = false; //cycle through electron filling for samples
		boolean ivCurve = false; //generate a curve of current versus voltage
		boolean multipleTemperatures = false; //cycle through multiple temperatures
		boolean multipleDisorders = true;
		//if none of the above are "true", then will simply go through single batch of samples with no frills
		
		if(multipleDiameters) {
			//The electron populations which will live on each nanoparticle. These correspond to a fixed volumetric density of .0015e-/nm^3 density with .5 ll
			List<Integer> nelecList = Arrays.asList(76, 109, 150, 199, 259, 329, 411, 506, 614, 736, 874, 1028, 1200, 1389, 1597, 1825, 2073);

			//The option exists to screen the charging energy
			for(double screeningFactor = 1.0; screeningFactor <= 1.0; screeningFactor += 1.0) {
				int j = 0;
				
				String title = "MIT Paper" + "\\" + Configuration.latticeStructure + diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_disorderType" + "percent2" +".txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("Diameter(nm) Mobility(cm^2/Vs) Pair_STD");
				
				for(double diam = 3; diam <= 8.5; diam = diam + 0.5) {
					System.out.println("nanoparticle diameter: " + String.valueOf(diam));
					double nelec_doub = nelecList.get(j)*.0016/.0015;
					int nelec = (int) Math.round(nelec_doub);
					System.out.println("Number of electrons is " + nelec);
					ExecutorService executor = Executors.newFixedThreadPool(4);
					List<Future<Double[]>> futures = new ArrayList<>();
					
					String folderName = numberOfNanoparticles +"_"+diam + "nm";
					if(preFlippedSamples) {
						for(int i=0; i<numberOfSamples; i++){
							futures.add(executor.submit(new Processor(folderName, i,false, T,nelec,screeningFactor, randomSeed, necking)));
						}
					}
					else {
						for(int i=0; i<numberOfSamples/2; i++){
							futures.add(executor.submit(new Processor(folderName, i,2*i, false, T,nelec,screeningFactor, randomSeed, necking)));
							futures.add(executor.submit(new Processor(folderName, i,2*i + 1, true, T,nelec,screeningFactor, randomSeed, necking)));
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
		
					//Process the results by averaging over the pairs of regular and inverted samples
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
					System.out.println(resultRaw);
					double[] resultProcessed = Utility.processResult(resultRaw);
					double total = 0.0;
					for(int i = 0; i < resultProcessed.length; i++) {
						total += resultProcessed[i];
					}
					
					//Get the mean, variance, and std
					double mean = total/resultProcessed.length;
					double variance = 0.0;
					for(int i=0; i< resultProcessed.length;i++) {
						variance += Math.pow((resultProcessed[i]-mean),2);
					}
					double std = Math.sqrt(variance/(resultProcessed.length-1));
					double std_mean =  std/Math.sqrt(resultProcessed.length);
	
					//Output the result
					String toAdd = String.valueOf(diam) + " " + mean + " " + std_mean;
					writer.println(toAdd);
					j += 1;
				}
				writer.close();
			}
		}
		
		else if(multipleDisorders) {
			// We have 450 (15x15x2) NPs, choose electron filling to be 0.5 e-/NP
			int nelec = 225;
			//The option exists to screen the charging energy
			for(double screeningFactor = 1.0; screeningFactor <= 1.0; screeningFactor += 1.0) {
				int j = 0;
				
				String title = "BinaryNSLpaper" + "\\" + Configuration.latticeStructure + diameter + "nm" +npArraySize + "_" + T + "K.txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("PercentDisorder Mobility(cm^2/Vs) Pair_STD");
				
				for(double pdisorder = .01; pdisorder <= .25; pdisorder = pdisorder + 0.01) {
					//double diam = 3.5; //nm
					System.out.println("nanoparticle % disorder: " + String.valueOf(pdisorder));
					System.out.println("Number of electrons is " + nelec);
					ExecutorService executor = Executors.newFixedThreadPool(16);
					List<Future<Double[]>> futures = new ArrayList<>();
					
					String folderName = numberOfNanoparticles +"_" + diameter + "nm_" +String.format("%.2f",pdisorder) + "D";
					if(preFlippedSamples) {
						for(int i=0; i<numberOfSamples; i++){
							futures.add(executor.submit(new Processor(folderName, i,false, T,nelec,screeningFactor, randomSeed, necking)));
						}
					}
					else {
						for(int i=0; i<numberOfSamples/2; i++){
							futures.add(executor.submit(new Processor(folderName, i,2*i, false, T,nelec,screeningFactor, randomSeed, necking)));
							futures.add(executor.submit(new Processor(folderName, i,2*i + 1, true, T,nelec,screeningFactor, randomSeed, necking)));
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
		
					//Process the results by averaging over the pairs of regular and inverted samples
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
					System.out.println(resultRaw);
					double[] resultProcessed = Utility.processResult(resultRaw);
					double total = 0.0;
					for(int i = 0; i < resultProcessed.length; i++) {
						total += resultProcessed[i];
					}
					
					//Get the mean, variance, and std
					double mean = total/resultProcessed.length;
					double variance = 0.0;
					for(int i=0; i< resultProcessed.length;i++) {
						variance += Math.pow((resultProcessed[i]-mean),2);
					}
					double std = Math.sqrt(variance/(resultProcessed.length-1));
					double std_mean =  std/Math.sqrt(resultProcessed.length);
	
					//Output the result
					String toAdd = String.valueOf(pdisorder) + " " + mean + " " + std_mean;
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
			
			//This is the on-site average electron filling
			for(double fraction = 0.8; fraction <= 1.0; fraction = fraction + 0.1) {
				int nelec = (int) Math.round(fraction*800);
				double diam = 6.5;
				System.out.println("Number Electrons: " + String.valueOf(nelec));
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
	
				//Run the simulations
				String folderName = numberOfNanoparticles +"_"+diam + "nm";
				if(preFlippedSamples) {
					for(int i=0; i<numberOfSamples; i++){
						futures.add(executor.submit(new Processor(folderName, i,false, T,nelec,screeningFactor, randomSeed, necking)));
					}
				}
				else {
					for(int i=0; i<numberOfSamples/2; i++){
						futures.add(executor.submit(new Processor(folderName, i,2*i, false, T,nelec,screeningFactor, randomSeed, necking)));
						futures.add(executor.submit(new Processor(folderName, i,2*i + 1, true, T,nelec,screeningFactor, randomSeed, necking)));
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
	
				//Process the results by averaging over the pairs of regular and inverted samples
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
				//System.out.println(resultRaw);
				double[] resultProcessed = Utility.processResult(resultRaw);
				
				//Calculate the mean mobility
				double total = 0.0;
				for(int i = 0; i < resultProcessed.length; i++) {
					total += resultProcessed[i];
				}
				double average = total/resultProcessed.length;
				
				//Output the result
				String toAdd = String.valueOf(nelec) + " " + average;
				writer.println(toAdd);
			}
			writer.close();
		}
		
		//This code is used for calculating the mobility gap for transport.
		else if(multipleTemperatures) {
			double screeningFactor = 1.0;
			
			//Simulate results for multiple disorder values
			List<String> disorderList = Arrays.asList(".00D",".015D", ".025D");
			for(String disorder: disorderList) {
				//String disorder = ".01D";
				String title = "TriDENS Results" + "/" + Configuration.latticeStructure + "_" + diameter + "nm" +npArraySize + "_" + eDensity + "_"  + "disorder" + disorder +".txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("Temperature(K) Mobility(cm^2/Vs) STD");
				
				//Necessary to look at multiple temperatures to calculate mobility gap
				for(double temp = 65.0; temp <= 80.0; temp += 2.5) {
					int nelec = 801; //1+delta e- per NP
					int crack = 0;
					double diam = 6.5;
					System.out.println("Temperature: " + temp);
					String folderName = numberOfNanoparticles +"_"+diam + "nm_200N_" + disorder;
					
					//Run the simulations
					ExecutorService executor = Executors.newFixedThreadPool(4);
					List<Future<Double[]>> futures = new ArrayList<>();
					if(preFlippedSamples) {
						for(int i=0; i<numberOfSamples; i++){
							futures.add(executor.submit(new Processor(folderName, i,false, temp,nelec,screeningFactor, randomSeed, necking)));
						}
					}
					else {
						for(int i=0; i<numberOfSamples/2; i++){
							futures.add(executor.submit(new Processor(folderName, i,2*i, false, temp,nelec,screeningFactor, randomSeed, necking)));
							futures.add(executor.submit(new Processor(folderName, i,2*i + 1, true, temp,nelec,screeningFactor, randomSeed, necking)));
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
		
					//Process the result, and average over regular and inverted samples
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
					}
					
					//Calculate the mean, variance, and standard deviation
					double average = total/nonZeroResults;
					double variance = 0.0;
					for(int i=0; i< nonZeroResults;i++) {
						variance += Math.pow((resultProcessed[i]-average),2);
					}
					double std = Math.sqrt(variance/(nonZeroResults-1));
					double std_mean =  std/Math.sqrt(nonZeroResults);
					
					//Print the results
					String toAdd = String.valueOf(temp) + " " + average +" " + std_mean;
					writer.println(toAdd);
				}
				writer.close();
			}
		}
		
		//This code is used for outputting the IV curve, useful for determining if we are in the linear regime
		else if(ivCurve) {
			double screeningFactor = 1.0;
			
			for(double voltageRatio = 0.25; voltageRatio <= 1.5; voltageRatio = voltageRatio + 0.25) {
				String title = "TriDENS Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_voltageRatio" + voltageRatio +".txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("SampleNumber Mobility(cm^2/Vs)");
				int nelec = 350;
				double diam = 6.5;
				String folderName = numberOfNanoparticles +"_"+diam + "nm";
				
				//Run the simulations
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
				if(preFlippedSamples) {
					for(int i=0; i<numberOfSamples; i++){
						futures.add(executor.submit(new Processor(folderName, i,i,false, T,nelec,screeningFactor, randomSeed, necking, voltageRatio)));
					}
				}
				else {
					for(int i=0; i<numberOfSamples/2; i++){
						futures.add(executor.submit(new Processor(folderName, i,2*i, false, T,nelec,screeningFactor, randomSeed, necking, voltageRatio)));
						futures.add(executor.submit(new Processor(folderName, i,2*i + 1, true, T,nelec,screeningFactor, randomSeed, necking, voltageRatio)));
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
	
				//Process the result and average over regular and inverted samples
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
				System.out.println(resultRaw);
				double[] resultProcessed = Utility.processResult(resultRaw);

				//Output the result
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
			String title = "TriDENS Results" + "\\" + Configuration.latticeStructure + "_" + diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K" + defect +".txt";
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs) Pair_STD");
			
			for(int k = 0; k < 1; k++) {
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
				String folderName = numberOfNanoparticles +"_"+diam + "nm";
				
				//submit jobs with an ordering id separate from their actual sample id
				//this will allow us to average over groups of runs on the same sample
				if(preFlippedSamples) {
					for(int i=0; i<numberOfSamples; i+=2){
						for(int j=0; j<repetitions; j++) {
							int nonflippedID = i*repetitions + 2*j;
							int flippedID = i*repetitions + 2*j + 1;

							futures.add(executor.submit(new Processor(folderName, i,nonflippedID, false, T,nelec,screeningFactor, randomSeed, necking)));
							futures.add(executor.submit(new Processor(folderName, i+1,flippedID, false, T,nelec, screeningFactor, randomSeed, necking)));
						}
					}
				}
				else {
					for(int i=0; i<numberOfSamples/2; i++){
						for(int j=0; j<repetitions; j++) {
							int nonflippedID = (2*i)*repetitions + 2*j;
							int flippedID = (2*i+1)*repetitions + 2*j + 1;
							
							futures.add(executor.submit(new Processor(folderName, 2*i,nonflippedID, false, T,nelec,screeningFactor, randomSeed, necking)));
							futures.add(executor.submit(new Processor(folderName, 2*i+1,flippedID, true, T,nelec, screeningFactor, randomSeed, necking)));
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
		
		//The default, most basic calculation. Simply calculate the mobility for pairs of samples. 
		else {
			int nelec = 800; //.002 e-/nm^3 triclinic 20x2x20
			int nhole = 0;
			double diam = 6.5;
			double screeningFactor = 1.0;
			String folderName = "20x4x20";
			String title = "New NP necking project" + "\\" + Configuration.latticeStructure + "_" + diameter + "nm_" +npArraySize + "_" + nelec + "e_" + T + "K_" + numberOfSamples+ "samples.txt";
			
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs)");
			
			//Run the simulations
			ExecutorService executor = Executors.newFixedThreadPool(4);
			List<Future<Double[]>> futures = new ArrayList<>();
			if(preFlippedSamples) {
				for(int i=0; i<numberOfSamples; i++){
					futures.add(executor.submit(new Processor(folderName, i,false, T,nelec,screeningFactor, randomSeed, necking)));
				}
			}
			else {
				for(int i=0; i<numberOfSamples/2; i++){
					futures.add(executor.submit(new Processor(folderName, i+sampleStart,2*i, false, T,nelec,screeningFactor, randomSeed, necking)));
					futures.add(executor.submit(new Processor(folderName, i+sampleStart,2*i + 1, true, T,nelec,screeningFactor, randomSeed, necking)));
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

			//Process the result and average over regular and inverted samples
			for(Future<Double[]> future: futures){
				try {
					resultRaw.add(future.get());
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}
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



