// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 16st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.numerics.SimpleVector.VectorNormType;

public class DarkFieldEllipsoid {

	
private SimpleVector scatterCoef;
private SimpleMatrix eigenVectors;
private SimpleMatrix scatterDirections;
private SimpleVector halfAxisLengths;
	

	/**
	 * This definition are taken from 2015-Vogel - X-Ray tensor tomography reconstruction
	 * paragraph 2.3 Constraint Enforcement and before
	 * @param scatterCoef - Reconstructed scattering coefficients (the n, not the ones
	 * where the square root was already taken)
	 * @param eigenValues - Eigenvalues of Ellipsoid
	 * @param eigenVectors - Eigenvector/ PCA Components of Ellipsoid
	 */
	public DarkFieldEllipsoid(SimpleMatrix scatterDirections, SimpleVector scatterCoef, SimpleVector eigenValues, SimpleMatrix eigenVectors){
		
		// Given by default over constructed
		this.scatterCoef = scatterCoef;
		this.eigenVectors = eigenVectors;
		this.scatterDirections = scatterDirections;
		
		// 
		
		// Average squared scattering magnitude
		double avgSqrScatMag = scatterCoef.norm(VectorNormType.VEC_NORM_L1)/scatterCoef.getLen();
		
		// Average eigenvalue magnitude
		double avgEigenValMag = eigenValues.norm(VectorNormType.VEC_NORM_L1)/eigenValues.getLen();
		
		// Define a size correction factor sigma for scaling the statistically defined ellipsoid to the point set S.
		double corrFactor = avgSqrScatMag*avgEigenValMag;
		
		// Allocate half-axis-length vector
		halfAxisLengths = new SimpleVector(3);
		
		/* 
		 * Calculating Ellipsoids Half-Axis-Length with respect to the orthonormal
		 * basic formed by the eigenVectors v1,v2,v3
		*/
		for(int i = 0; i < 3; i++){
			double val = Math.sqrt(corrFactor*eigenValues.getElement(i));
			halfAxisLengths.setElementValue(i, val);
		}
		
	}
	
	
	
	
	/**
	 * Transforms one point into the ellipsoid's coordinate frame
	 * thus obtaining a vector [x,y,z] = V'*r
	 * @param myPoint - has to be of dimension 3
	 * @return
	 */
	private SimpleVector transformPointIntoEigenVectorSystem(SimpleVector myPoint){
		SimpleVector transformed = transformPointIntoEigenVectorSystem(myPoint, eigenVectors);
		return transformed;
	}
	
	/**
	 * Transforms one point into the ellipsoid's coordinate frame
	 * thus obtaining a vector [x,y,z] = V'*r
	 * @param myPoint - has to be of dimension 3
	 * @return
	 */
	public static SimpleVector transformPointIntoEigenVectorSystem(SimpleVector myPoint, SimpleMatrix eigenVectors){
		SimpleVector transformed = SimpleOperators.multiply(eigenVectors.transposed(), myPoint);
		return transformed;
	}
	
	
	/**
	 * Consider that squared coefficients are reconstructed, the correctly projected coefficients are
	 * given by n = sigma^2
	 * @return
	 */
	public SimpleVector calculateSquaredProjectedCoefficients(){

		/* Calculate the contribution (projection) of every
		* scatter direction onto the ellipsoid 
		*/
		
		SimpleVector projScatterCoef = new SimpleVector(this.scatterCoef.getLen());
		
		for(int channel = 0; channel < scatterCoef.getLen(); channel++){
		
		// Get the scatter direction of the current scatter vector.
		// Remember, scatter directions are saved in the columns
		SimpleVector curDir = scatterDirections.getCol(channel);
			
		// Calculate vector in ellipsoid coordinate frame
		SimpleVector point = transformPointIntoEigenVectorSystem(curDir);
		
		double curCoef;
		// Scale points by the half axis lengths like eq. 18 in vogel's paper
		if(halfAxisLengths.min()==0){
			curCoef = 0;
		} else {
		SimpleVector scaledPoint = SimpleOperators.divideElementWise(point,halfAxisLengths);
		 curCoef = 1/(scaledPoint.normL2()*scaledPoint.normL2());
		}
		
		projScatterCoef.setElementValue(channel, curCoef);
		}
		
		return projScatterCoef;
	}
	
	
	
}
	
