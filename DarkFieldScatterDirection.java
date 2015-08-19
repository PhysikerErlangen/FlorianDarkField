package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldScatterDirection {

	public static SimpleMatrix getScatterDirectionMatrix(int numScatterVectors){

		// Create matrix, 3 rows because of 3 world coordinates x,y,z
		// and numScatterVectors columns
		SimpleMatrix scatterMatrix = new SimpleMatrix(3,numScatterVectors);
		
		// For the first test implementations
	if(numScatterVectors == 3){
		
		scatterMatrix.setColValue(0, new SimpleVector(1f,0f,0f));
		scatterMatrix.setColValue(1, new SimpleVector(0f,1f,0f));
		scatterMatrix.setColValue(2, new SimpleVector(0f,0f,1f));
	}	
	
	if(numScatterVectors == 13){
		
		scatterMatrix.setColValue(0, new SimpleVector(1f,0f,0f).normalizedL2());
		scatterMatrix.setColValue(1, new SimpleVector(0f,1f,0f).normalizedL2());
		scatterMatrix.setColValue(2, new SimpleVector(0f,0f,1f).normalizedL2());
		scatterMatrix.setColValue(3, new SimpleVector(1f,1f,1f).normalizedL2());
		scatterMatrix.setColValue(4, new SimpleVector(1f,1f,-1f).normalizedL2());
		scatterMatrix.setColValue(5, new SimpleVector(1f,-1f,1f).normalizedL2());
		scatterMatrix.setColValue(6, new SimpleVector(-1f,1f,1f).normalizedL2());
		
		scatterMatrix.setColValue(7, new SimpleVector(1f,1f,0f).normalizedL2());
		scatterMatrix.setColValue(8, new SimpleVector(-1f,1f,0f).normalizedL2());
		
		scatterMatrix.setColValue(9, new SimpleVector(0f,1f,1f).normalizedL2());
		scatterMatrix.setColValue(10, new SimpleVector(0f,-1f,1f).normalizedL2());
		
		scatterMatrix.setColValue(11, new SimpleVector(1f,0f,1f).normalizedL2());
		scatterMatrix.setColValue(12, new SimpleVector(-1f,0f,1f).normalizedL2());
	}
	
	else if(numScatterVectors == 1){
		scatterMatrix.setColValue(0, new SimpleVector(1f));
	}
	
	return scatterMatrix;
	}
	
	
}
