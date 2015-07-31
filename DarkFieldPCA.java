package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.geometry.shapes.activeshapemodels.PCA;
import edu.stanford.rsl.conrad.geometry.shapes.mesh.DataMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldPCA{

	// Contains the different scatterDirections (is constant and given as input)
	ArrayList<SimpleVector> scatterDirections;
	
	// Scatter points contains the scatter directions of the respective 
	// voxel which are alread weighted by the earlier calculated
	// scatter coefficients. These points are then used for calculating
	// the principal axes by PCA
	ArrayList<SimpleVector> scatterPoints;
	
	PCA myPCA;
	
	public DarkFieldPCA (){
		
	
	}
	
	private void calcDataMatrix(){
		
		
		
	}
	
	// Calculates the different scatter directions with their corresponding weights. 
	private void calcSetOfScatterPoints(SimpleVector myScatterWeights){
		
		scatterPoints = new ArrayList<SimpleVector>();
		
		int numScatterVectors = myScatterWeights.getLen();
		
		for (int channel = 0; channel < numScatterVectors; channel++){
			
			double myWeight = myScatterWeights.getElement(channel);
			SimpleVector vec1 = scatterDirections.get(channel);
			vec1.multipliedBy(myWeight);
			SimpleVector vec2 = vec1.multipliedBy(-1f);
			
			scatterPoints.add(vec1);
			scatterPoints.add(vec2);
			
		}
		
		
	}
	
	
}
