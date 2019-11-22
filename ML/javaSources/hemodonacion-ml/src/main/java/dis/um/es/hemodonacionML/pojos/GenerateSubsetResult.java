package dis.um.es.hemodonacionML.pojos;

import java.io.Serializable;

import weka.attributeSelection.AttributeSelection;
import weka.core.Instances;

public class GenerateSubsetResult implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1242505190513835434L;
	private Instances reducedDataset;
	private String name;
	private AttributeSelection attributeSelector;
	public GenerateSubsetResult(Instances reducedDataset, String name, AttributeSelection attributeSelector) {
		super();
		this.reducedDataset = reducedDataset;
		this.name = name;
		this.attributeSelector = attributeSelector;
	}
	public Instances getReducedDataset() {
		return reducedDataset;
	}
	public void setReducedDataset(Instances reducedDataset) {
		this.reducedDataset = reducedDataset;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public AttributeSelection getAttributeSelector() {
		return attributeSelector;
	}
	public void setAttributeSelector(AttributeSelection attributeSelector) {
		this.attributeSelector = attributeSelector;
	}
	
	
}
