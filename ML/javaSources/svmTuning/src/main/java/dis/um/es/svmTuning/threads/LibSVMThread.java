package dis.um.es.svmTuning.threads;

import java.util.Random;
import java.util.concurrent.Callable;

import dis.um.es.svmTuning.pojos.TestClassifierResult;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.LibSVM;
import weka.core.Instances;

public class LibSVMThread implements Callable<TestClassifierResult>{
	private Instances dataset;
	private LibSVM libSVM;
	private double gamma;
	private double cost;
	private int iteration;
	private int crossValidationFolds;
	public LibSVMThread(Instances dataset, double gamma, double cost, int iteration, int crossValidationFolds) {
		super();
		this.gamma = gamma;
		this.cost = cost;
		this.iteration = iteration;
		this.dataset = new Instances(dataset);
		this.crossValidationFolds = crossValidationFolds;
		libSVM = new LibSVM();
		libSVM.setCost(cost);
		libSVM.setGamma(gamma);
		libSVM.setSeed(iteration);
	}
	public LibSVM getLibSVM() {
		return libSVM;
	}
	public void setLibSVM(LibSVM libSVM) {
		this.libSVM = libSVM;
	}
	public double getGamma() {
		return gamma;
	}
	public void setGamma(double gamma) {
		this.gamma = gamma;
	}
	public double getCost() {
		return cost;
	}
	public void setCost(double cost) {
		this.cost = cost;
	}
	public int getIteration() {
		return iteration;
	}
	public void setIteration(int iteration) {
		this.iteration = iteration;
	}
	public Instances getDataset() {
		return dataset;
	}
	public void setDataset(Instances dataset) {
		this.dataset = dataset;
	}
	public TestClassifierResult call() throws Exception {
		Evaluation evaluation = new Evaluation(dataset);
		evaluation.crossValidateModel(libSVM, dataset, crossValidationFolds, new Random(iteration));
		return new TestClassifierResult(evaluation, libSVM, dataset, iteration);
	}
	
	
}
