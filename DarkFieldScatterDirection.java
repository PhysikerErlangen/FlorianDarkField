package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldScatterDirection {

	public static SimpleMatrix getScatterDirectionMatrix(int numScatterVectors){
		return calcScatterMatrixBySpiralMethod(numScatterVectors);
	}
	
	
	public static SimpleMatrix predefinedScatterMatrix(int numScatterVectors){

		// Create matrix, 3 rows because of 3 world coordinates x,y,z
		// and numScatterVectors columns
		SimpleMatrix scatterMatrix = new SimpleMatrix(3,numScatterVectors);
		
		// For the first test implementations
	if(numScatterVectors == 3){
		
		scatterMatrix.setColValue(0, new SimpleVector(1f,0f,0f));
		scatterMatrix.setColValue(1, new SimpleVector(0f,1f,0f));
		scatterMatrix.setColValue(2, new SimpleVector(0f,0f,1f));
	}	
	else if(numScatterVectors == 13){
		
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
else if(numScatterVectors == 7){
		
		scatterMatrix.setColValue(0, new SimpleVector(1f,0f,0f).normalizedL2());
		scatterMatrix.setColValue(1, new SimpleVector(0f,1f,0f).normalizedL2());
		scatterMatrix.setColValue(2, new SimpleVector(0f,0f,1f).normalizedL2());
		scatterMatrix.setColValue(3, new SimpleVector(1f,1f,1f).normalizedL2());
		scatterMatrix.setColValue(4, new SimpleVector(1f,1f,-1f).normalizedL2());
		scatterMatrix.setColValue(5, new SimpleVector(1f,-1f,1f).normalizedL2());
		scatterMatrix.setColValue(6, new SimpleVector(-1f,1f,1f).normalizedL2());
	}
	
	else if(numScatterVectors == 1){
		scatterMatrix.setColValue(0, new SimpleVector(1f));
	}
	
	return scatterMatrix;
	}
	
	
	
	/**
	 * More information can be found here:
	 * http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere
	 * and the following paper:
	 *  Distributing many points on a sphere
	 *  by E. B. Saff, A. B. J. Kuijlaars
	 * LINK to Paper: http://link.springer.com/article/10.1007%2FBF03024331#page-1
	 * @param numScatterVectors
	 * @return
	 */
	private static SimpleMatrix calcScatterMatrixBySpiralMethod(int numScatterVectors){
		
		// Create matrix, 3 rows because of 3 world coordinates x,y,z
		// and numScatterVectors columns
		SimpleMatrix scatterMatrix = new SimpleMatrix(3,numScatterVectors);
		
		double s = 3.6/Math.sqrt(numScatterVectors);
		/*
		 * Caution: Usually 2.0/N, but we only want to sample one half of the sphere
		 * so take half the value
		 */
		double dz = 1.0/numScatterVectors;
		
		double l = 0;
		
		double z = 1 - dz/2.0;
		
		for(int k = 0; k < numScatterVectors; k++){
			
			double r = Math.sqrt(1-z*z);
			SimpleVector vec = new SimpleVector(Math.cos(l)*r,Math.sin(l)*r,z);
			z = z - dz;
			l = l + s/r;
			
			
			scatterMatrix.setColValue(k, vec);
			
		
		}
		
		return scatterMatrix;
		
	}
	
}
