// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 16st, 2015
//

package edu.stanford.rsl.science.darkfield.Florian;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

import edu.stanford.rsl.conrad.numerics.Solvers;

public class DarkFieldGaussian {

	// Covariance matrix
	SimpleMatrix cov;
	
	// Mean vector
	SimpleVector mu;
	
	//TODO Implement this for N-dimension
	public DarkFieldGaussian(double sigma1, double sigma2, double sigma3){
		SimpleVector diag = new SimpleVector(sigma1*sigma1,sigma2*sigma2,sigma3*sigma3);		
		init(diag);
		
	}
	
	public DarkFieldGaussian(SimpleVector diag){
		init(diag);
	}
	
	private void init(SimpleVector diag){
		cov = new SimpleMatrix(3,3);
		cov.setDiagValue(diag);
	}
	
	
	public double calcGaussian(SimpleVector x){
		
		SimpleVector help = Solvers.solveLinearSysytemOfEquations(cov, x);
		//SimpleVector help = SimpleOperators.multiply(cov.inverse(INVERT_SVD), x);
		double exponent = SimpleOperators.multiplyInnerProd(x, help);
		double gauss = Math.exp(-0.5*exponent);
		
		return gauss;
		
	}
	
	
}
