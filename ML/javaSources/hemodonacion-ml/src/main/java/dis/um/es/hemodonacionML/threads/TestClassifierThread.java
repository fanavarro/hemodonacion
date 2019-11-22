package dis.um.es.hemodonacionML.threads;


import java.util.Random;
import java.util.concurrent.Callable;

import dis.um.es.hemodonacionML.pojos.TestClassifierResult;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.core.Instances;

public class TestClassifierThread implements Callable<TestClassifierResult>{
	private Classifier classifier;
	private Instances dataset;
	private String datasetName;
	private int crossValidationFolds;
	private int iteration;
	public TestClassifierThread(Classifier classifier, Instances dataset, String datasetName, int crossValidationFolds, int iteration) throws Exception {
		super();
		this.classifier = AbstractClassifier.makeCopy(classifier);
		this.dataset = new Instances(dataset);
		this.crossValidationFolds = crossValidationFolds;
		this.datasetName = datasetName;
		this.iteration = iteration;
	}
	public Classifier getClassifier() {
		return classifier;
	}
	public void setClassifier(Classifier classifier) {
		this.classifier = classifier;
	}
	public Instances getDataset() {
		return dataset;
	}
	public void setDataset(Instances dataset) {
		this.dataset = dataset;
	}
	public int getCrossValidationFolds() {
		return crossValidationFolds;
	}
	public void setCrossValidationFolds(int crossValidationFolds) {
		this.crossValidationFolds = crossValidationFolds;
	}
	
	
	public String getDatasetName() {
		return datasetName;
	}
	public void setDatasetName(String datasetName) {
		this.datasetName = datasetName;
	}
	@Override
	public TestClassifierResult call() throws Exception {
		Evaluation evaluation = new Evaluation(dataset);
		evaluation.crossValidateModel(classifier, dataset, crossValidationFolds, new Random(iteration));
		return new TestClassifierResult(evaluation, datasetName, classifier, dataset, iteration);
	}
	
	
}
