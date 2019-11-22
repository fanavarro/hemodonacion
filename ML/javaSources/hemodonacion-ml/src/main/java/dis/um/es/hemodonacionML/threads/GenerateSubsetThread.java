package dis.um.es.hemodonacionML.threads;

import java.util.concurrent.Callable;

import dis.um.es.hemodonacionML.pojos.GenerateSubsetResult;
import weka.attributeSelection.ASEvaluation;
import weka.attributeSelection.ASSearch;
import weka.attributeSelection.AttributeSelection;
import weka.attributeSelection.WrapperSubsetEval;
import weka.core.Instances;

public class GenerateSubsetThread implements Callable<GenerateSubsetResult>{
	private int SEED = 12345;
	private ASEvaluation evaluationMethod;
	private ASSearch searchMethod;
	private Instances data;
	
	
	
	public GenerateSubsetThread(ASEvaluation evaluationMethod, ASSearch searchMethod, Instances data) throws Exception {
		super();
		this.evaluationMethod = ASEvaluation.makeCopies(evaluationMethod, 1)[0];
		this.searchMethod = ASSearch.makeCopies(searchMethod, 1)[0];
		this.data = new Instances(data);
	}


	private static String getName(ASEvaluation evaluation, ASSearch search){
		String name = evaluation.getClass().getSimpleName();
		if(evaluation instanceof WrapperSubsetEval){
			name += String.format("(%s)", ((WrapperSubsetEval)evaluation).getClassifier().getClass().getSimpleName()); 
		}
		name += "-" + search.getClass().getSimpleName();
		return name;
	}



	@Override
	public GenerateSubsetResult call() throws Exception {
		String name = getName(evaluationMethod, searchMethod);
		AttributeSelection attributeSelector = new AttributeSelection();
		attributeSelector.setSeed(SEED);
		attributeSelector.setEvaluator(evaluationMethod);
		attributeSelector.setSearch(searchMethod);
		attributeSelector.SelectAttributes(data);
		Instances reducedInstances = attributeSelector.reduceDimensionality(data);
		return new GenerateSubsetResult(reducedInstances, name, attributeSelector);
	}
}
