package dis.um.es.svmTuning.pojos;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import weka.classifiers.Evaluation;
import weka.classifiers.functions.LibSVM;
import weka.core.Instances;

public class TestClassifierResult {
	private Evaluation evaluation;
	private LibSVM classifier;
	private Instances dataset;
	private int iteration;

	public TestClassifierResult(Evaluation evaluation, LibSVM classifier, Instances dataset, int iteration) {
		super();
		this.evaluation = evaluation;
		this.classifier = classifier;
		this.dataset = dataset;
		this.iteration = iteration;
	}

	public Evaluation getEvaluation() {
		return evaluation;
	}

	public void setEvaluation(Evaluation evaluation) {
		this.evaluation = evaluation;
	}


	public LibSVM getClassifier() {
		return classifier;
	}

	public void setClassifier(LibSVM classifier) {
		this.classifier = classifier;
	}

	public Instances getDataset() {
		return dataset;
	}

	public void setDataset(Instances dataset) {
		this.dataset = dataset;
	}

	@Override
	public String toString() {
		DecimalFormatSymbols symbols = new DecimalFormatSymbols(new Locale("en", "UK"));
		DecimalFormat df = new DecimalFormat("0.000", symbols);
		StringBuilder sb = new StringBuilder();

		sb.append(classifier.getClass().getSimpleName()).append("\t");
		
		sb.append(df.format(classifier.getGamma())).append("\t");
		
		sb.append(df.format(classifier.getCost())).append("\t");
		
		sb.append(iteration).append("\t");

		sb.append(df.format(evaluation.correct())).append("\t");

		sb.append(df.format(evaluation.pctCorrect())).append("\t");

		sb.append(df.format(evaluation.incorrect())).append("\t");

		sb.append(df.format(evaluation.pctIncorrect())).append("\t");

		sb.append(df.format(evaluation.kappa())).append("\t");

		sb.append(df.format(evaluation.errorRate())).append("\t");

		for (int i = 0; i < dataset.attribute(dataset.classIndex()).numValues(); i++) {
			sb.append(df.format(evaluation.numTruePositives(i))).append("\t");
			sb.append(df.format(evaluation.truePositiveRate(i))).append("\t");
			sb.append(df.format(evaluation.numTrueNegatives(i))).append("\t");
			sb.append(df.format(evaluation.trueNegativeRate(i))).append("\t");
			sb.append(df.format(evaluation.numFalsePositives(i))).append("\t");
			sb.append(df.format(evaluation.falsePositiveRate(i))).append("\t");
			sb.append(df.format(evaluation.numFalseNegatives(i))).append("\t");
			sb.append(df.format(evaluation.falseNegativeRate(i))).append("\t");
			sb.append(df.format(evaluation.recall(i))).append("\t");
			sb.append(df.format(evaluation.precision(i))).append("\t");
			sb.append(df.format(evaluation.fMeasure(i))).append("\t");
			sb.append(df.format(evaluation.areaUnderROC(i))).append("\t");
		}
		sb.append(df.format(evaluation.numInstances())).append("\n");

		return sb.toString();
	}

	public String getHeader() {
		StringBuilder header = new StringBuilder();
		header.append("Classifier\t");
		header.append("Gamma\t");
		header.append("Cost\t");
		header.append("Iteration\t");
		header.append("Correct\t");
		header.append("Percentage Correct\t");
		header.append("Incorrect\t");
		header.append("Percentage Incorrect\t");
		header.append("Kappa\t");
		header.append("Error Rate\t");
		for (int i = 0; i < dataset.attribute(dataset.classIndex()).numValues(); i++) {
			String classValue = String.format("(%s)\t", dataset.attribute(dataset.classIndex()).value(i));
			header.append("True Positives ").append(classValue);
			header.append("True Positive Rate ").append(classValue);
			header.append("True Negatives ").append(classValue);
			header.append("True Negative Rate ").append(classValue);
			header.append("False Positives ").append(classValue);
			header.append("False Positive Rate ").append(classValue);
			header.append("False Negatives ").append(classValue);
			header.append("False Negative Rate ").append(classValue);
			header.append("Recall ").append(classValue);
			header.append("Precision ").append(classValue);
			header.append("F-Measure ").append(classValue);
			header.append("ROC Area ").append(classValue);
		}
		header.append("Total Instancies\n");
		return header.toString();
	}

}
