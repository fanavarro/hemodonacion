package dis.um.es.svmTuning.main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import dis.um.es.svmTuning.pojos.TestClassifierResult;
import dis.um.es.svmTuning.threads.LibSVMThread;
import dis.um.es.svmTuning.utils.Utils;
import weka.core.Instances;

public class Main {

	@Option(name = "-f", required=true, usage = "ARFF file to evaluate SVM classifiers.", metaVar="INPUT")
	private File inputFile;
	
	@Option(name = "-t", required=false, usage = "Number of cores to use.")
	private int nCores = 1;
	
	@Option(name = "-i", required=false, usage = "Iterations for each classifier.")
	private int iterations = 1;
	
	private int gammaLowerBound = -15;
	private int gammaUpperBound = 3;
	private int costLowerBound = -5;
	private int costUpperBound = 15;
	
	private int crossValidationFolds = 10;
	
	
	public static void main(String[] args) throws Exception {
		new Main().doMain(args);
	}
	
	public void doMain(String[] args) throws CmdLineException, IOException, InterruptedException, ExecutionException {
		CmdLineParser parser = new CmdLineParser(this);
		
		parser.parseArgument(args);
		
		BufferedReader dataFile = Utils.readDataFile(inputFile);
		Instances data = new Instances(dataFile);
		data.setClassIndex(data.numAttributes() - 1);

		
		
		/* Preparar y lanzar hilos */
		List<LibSVMThread> threads = this.generateThreads(data);
		System.out.println("Lanzando trabajo");
		String csv = runThreads(threads);
		
		writeToFile(inputFile.getParentFile().getAbsolutePath() + "/SVMEvaluationResults.tsv", csv);
		System.out.println("FIN");
	}
	
	private void writeToFile(String filename, String s) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
		writer.write(s);
		writer.close();
	}
	
	private String runThreads(List<LibSVMThread> threads) throws InterruptedException, ExecutionException{
		StringBuilder sb = new StringBuilder();
		ExecutorService executor = Executors.newFixedThreadPool(this.nCores);
		List<Future<TestClassifierResult>> futures = new LinkedList<Future<TestClassifierResult>>();
		for(Callable<TestClassifierResult> thread : threads){
			futures.add(executor.submit(thread));
		}
		
		/* Recuperar resultados */
		boolean headerIncluded = false;
		for(Future<TestClassifierResult> future : futures){
			TestClassifierResult result = future.get();
			if(!headerIncluded){
				sb.append(result.getHeader());
				headerIncluded = true;
			}
			sb.append(result.toString());
		}
		executor.shutdown();
		
		return sb.toString();
	}

	private List<LibSVMThread> generateThreads(Instances dataset) {
		List<LibSVMThread> threads = new ArrayList<LibSVMThread> ();
		for(int gamma = gammaLowerBound ; gamma <= gammaUpperBound; gamma++){
			for(int cost = costLowerBound; cost <= costUpperBound; cost++){
				for(int iteration = 0; iteration < iterations; iteration++){
					threads.add(new LibSVMThread(dataset, Math.pow(2, gamma), Math.pow(2, cost), iteration, crossValidationFolds));
				}
			}
		}
		return threads;
	}

}
