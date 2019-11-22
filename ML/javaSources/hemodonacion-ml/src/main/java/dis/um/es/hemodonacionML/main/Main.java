package dis.um.es.hemodonacionML.main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import dis.um.es.hemodonacionML.pojos.GenerateSubsetResult;
import dis.um.es.hemodonacionML.pojos.TestClassifierResult;
import dis.um.es.hemodonacionML.threads.GenerateSubsetThread;
import dis.um.es.hemodonacionML.threads.TestClassifierThread;
import dis.um.es.hemodonacionML.utils.Utils;
import weka.attributeSelection.ASEvaluation;
import weka.attributeSelection.ASSearch;
import weka.attributeSelection.BestFirst;
import weka.attributeSelection.CfsSubsetEval;
import weka.attributeSelection.GreedyStepwise;
import weka.attributeSelection.MultiObjectiveEvolutionarySearch;
import weka.attributeSelection.PSOSearch;
import weka.attributeSelection.WrapperSubsetEval;
import weka.classifiers.Classifier;
import weka.classifiers.bayes.BayesNet;
import weka.classifiers.evaluation.EvaluationUtils;
import weka.classifiers.functions.LibSVM;
import weka.classifiers.functions.MultilayerPerceptron;
import weka.classifiers.lazy.IBk;
import weka.classifiers.lazy.KStar;
import weka.classifiers.rules.PART;
import weka.classifiers.trees.J48;
import weka.classifiers.trees.RandomForest;
import weka.core.Instances;

public class Main {
	
	private static Classifier[] CLASSIFIERS = { new BayesNet(), new J48(), new RandomForest(), new PART(), new IBk(),
			new KStar(), new LibSVM(), new MultilayerPerceptron() };
	
	private static ASEvaluation [] ATTRIBUTE_SELECTION_EVALUATORS = generateAttributeSelectionEvaluators();
	private static ASSearch [] ATTRIBUTE_SELECTION_SEARCH_METHOD = {new BestFirst(), new GreedyStepwise(), new MultiObjectiveEvolutionarySearch(), new PSOSearch()};



	@Option(name = "-f", required=true, usage = "Original ARFF file to extract subsets and evaluate classifiers.", metaVar="INPUT")
	private File inputFile;
	
	@Option(name = "-t", required=false, usage = "Number of cores to use.")
	private int nCores = 1;
	
	@Option(name = "-i", required=false, usage = "Iterations for each classifier.")
	private int iterations = 1;
	
	public static void main(String[] args) throws Exception {
		new Main().doMain(args);
	}
	
	public void doMain(String[] args) throws Exception {

		CmdLineParser parser = new CmdLineParser(this);
		
		parser.parseArgument(args);
		
		EvaluationUtils eu = new EvaluationUtils();
		eu.setSeed(12345);
		
		BufferedReader dataFile = Utils.readDataFile(inputFile);
		Instances data = new Instances(dataFile);
		data.setClassIndex(data.numAttributes() - 1);

		System.out.println("GENERATING SUBSETS");
		Map<String, Instances> subsets = generateSubsets(data);
		writeSubsetSummary(inputFile.getParentFile().getAbsolutePath() + "/subsetsSummary.txt", subsets);
		
		System.out.println("SAVING SUBSETS INTO DISK");
		writeArff(subsets, inputFile.getParentFile().getAbsolutePath() + "/subsets");
		
		System.out.println("TESTING CLASSIFIERS");
		subsets.put("original-original", data); // Evaluation of original dataset.
		String result = testClassifiers(Arrays.asList(CLASSIFIERS), subsets);
		writeToFile(inputFile.getParentFile().getAbsolutePath() + "/evaluationResults.tsv", result);
	}
	
	private void writeSubsetSummary(String outputFile, Map<String, Instances> subsets) throws IOException {
		StringBuilder sb = new StringBuilder();
		for(Entry <String, Instances> entry : subsets.entrySet()){
			sb.append(entry.getKey()).append('\t');
			for(int i = 0; i < entry.getValue().numAttributes(); i++){
				sb.append(entry.getValue().attribute(i).name()).append(',');
			}
			sb.append('\n');
		}
		writeToFile(outputFile, sb.toString());
	}

	private void writeArff(Map<String, Instances> subsets, String outputFolderPath) throws IOException{
		File outputFolder = new File(outputFolderPath);
		if(!outputFolder.exists() || !outputFolder.isDirectory()){
			outputFolder.mkdirs();
		}
		for(Entry<String, Instances> entry : subsets.entrySet()){
			String subsetFilename = entry.getKey() + ".arff";
			Instances subset = entry.getValue();
			writeToFile(outputFolder.getAbsolutePath() + "/" + subsetFilename, subset.toString());
		}
	}
	
	private void writeToFile(String filename, String s) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
		writer.write(s);
		writer.close();
	}

	private static ASEvaluation[] generateAttributeSelectionEvaluators(){
		List<ASEvaluation> evaluators = new ArrayList<ASEvaluation>();
		evaluators.add(new CfsSubsetEval());
		
		WrapperSubsetEval evaluator = new WrapperSubsetEval();
		evaluator.setClassifier(new LibSVM());
		evaluators.add(evaluator);
		
		evaluator = new WrapperSubsetEval();
		evaluator.setClassifier(new PART());
		evaluators.add(evaluator);
		
		evaluator = new WrapperSubsetEval();
		evaluator.setClassifier(new J48());
		evaluators.add(evaluator);
		
		evaluator = new WrapperSubsetEval();
		evaluator.setClassifier(new KStar());
		evaluators.add(evaluator);
		
		evaluator = new WrapperSubsetEval();
		evaluator.setClassifier(new BayesNet());
		evaluators.add(evaluator);
		
		ASEvaluation [] arrayEvaluations = new ASEvaluation[evaluators.size()];
		return evaluators.toArray(arrayEvaluations);
	}
	
	private Map<String, Instances> generateSubsets(Instances data) throws Exception{
		Map<String, Instances> subsets = new HashMap<String, Instances>();
//		ExecutorService executor = Executors.newSingleThreadExecutor();
		ExecutorService executor = Executors.newFixedThreadPool(this.nCores);
		List<Future<GenerateSubsetResult>> futures = new LinkedList<Future<GenerateSubsetResult>>();
		/* Preparar y lanzar los hilos */
		for(ASEvaluation evaluationMethod : ATTRIBUTE_SELECTION_EVALUATORS){
			for (ASSearch searchMethod : ATTRIBUTE_SELECTION_SEARCH_METHOD){
				GenerateSubsetThread thread = new GenerateSubsetThread(evaluationMethod, searchMethod, data);
				futures.add(executor.submit(thread));
			}
		}
		
		/* Recuperar los resultados */
		for(Future<GenerateSubsetResult> future : futures){
			GenerateSubsetResult result = future.get();
			subsets.put(result.getName(), result.getReducedDataset());
			System.out.println(result.getName() + "\n" + result.getAttributeSelector().toResultsString());
		}
		executor.shutdown();
		return subsets;
	}

	private String testClassifiers(List<Classifier> classifiers, Map<String, Instances> datasources) throws Exception{
		StringBuilder sb = new StringBuilder();
		int crossValidationFolds = 10;
		
		/* Preparar y lanzar hilos */
		ExecutorService executor = Executors.newFixedThreadPool(this.nCores);
		List<Future<TestClassifierResult>> futures = new LinkedList<Future<TestClassifierResult>>();
		for(Entry<String, Instances> datasource : datasources.entrySet()){
			Instances dataset = datasource.getValue();
			String datasetName = datasource.getKey();
			dataset.setClassIndex(dataset.numAttributes() - 1);

			for (Classifier classifier : CLASSIFIERS) {
				for(int i = 0; i < this.iterations; i++){
					TestClassifierThread thread = new TestClassifierThread(classifier, dataset, datasetName, crossValidationFolds, i);
					futures.add(executor.submit(thread));
				}
			}
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
}
