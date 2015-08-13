/**
 * 
 */
package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix.MatrixNormType;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.numerics.SimpleVector.VectorNormType;

/**
 * @author schiffers
 *
 */
public class DarkFieldPCA_Test {

	DarkFieldPCA testPCA;
	
	SimpleMatrix eigenVecFromMatlab;
	SimpleVector eigenValuesFromMatlab;

	/**
	 * @throws java.lang.Exception
	 */
	/**
	 * @throws Exception
	 */
	@Before
	public void setUp() throws Exception {
		
		// Values taken from example out of corresponding Matlab file
		double[][] data = new double[][]{
				   {-0.4615,  -22.1538,   -5.7692},
				   {-6.4615 , -19.1538  ,  3.2308},
				   {3.5385 ,   7.8462  , -3.7692},
				   {3.5385 , -17.1538  , -3.7692},
				   {-0.4615 ,   3.8462 ,  -5.7692},
				   {3.5385  ,  6.8462  , -2.7692},
				   {-4.4615 ,  22.8462 ,   5.2308},
				   {-6.4615 , -17.1538 ,  10.2308},
				   {-5.4615 ,   5.8462 ,   6.2308},
				   {13.5385 ,  -1.1538 ,  -7.7692},
				   {-6.4615 ,  -8.1538 ,  11.2308},
				   { 3.5385 ,  17.8462 ,  -2.7692},
				   {2.5385  , 19.8462 ,  -3.7692},
				  
				}; 
	
		// Data PC components should be (see matlab)
		double[][] eigenVecMatlab = new double[][]{
				{0.1105 ,  -0.6528 ,   0.7494},
			    {0.9903  ,  0.1360  , -0.0276},
			    {-0.0839 ,   0.7452 ,   0.6615}
		};
		eigenVecFromMatlab = new SimpleMatrix(eigenVecMatlab);
		// Eigenvalues calculated by matlab
		eigenValuesFromMatlab = new SimpleVector(  245.6525,  65.6918,  6.4250);
		
		int n = data.length;
		
		SimpleVector myScatterCoef = new SimpleVector(n);
		for(int i = 0; i < n; i++){
			myScatterCoef.setElementValue(i, 1);
		}
		
		SimpleMatrix scatterDirections = new SimpleMatrix(data);
		scatterDirections = scatterDirections.transposed();
	
		testPCA = new DarkFieldPCA(scatterDirections, myScatterCoef);
		
		testPCA.run();
		
	}



	@Test
	public void testEigenVectors() {
		SimpleMatrix eigenVectors = testPCA.getEigenVectors();
		
		System.out.println("EigenVectors are: "+eigenVectors.transposed());
		
		SimpleMatrix diffMatrix = SimpleOperators.subtract(eigenVectors, eigenVecFromMatlab);

		
		
		System.out.println("Difference between real Eigenvectors is: ");
		System.out.println(eigenVectors.transposed());
		
		double myNorm = diffMatrix.norm(MatrixNormType.MAT_NORM_FROBENIUS);
		System.out.println("The norm of differences is: " + myNorm);
		
		double th = 0.1;
		
		
		assertTrue("Norm of Eigenvector difference matrix is too large",myNorm < th);
	}
	
	@Test
	public void testEigenValues() {
		
		SimpleVector eigenValues = 	testPCA.getEigenValues();
		System.out.println("EigenValues in descending order: " + eigenValues);
		
		SimpleVector diffMatrix = SimpleOperators.subtract(eigenValues,eigenValuesFromMatlab);
		
		double norm = diffMatrix.norm(VectorNormType.VEC_NORM_L1);
		double th = 20;
		System.out.println("Norm of differnces: " +norm);
		assertTrue("This works fine",norm < th);
	}
	
	@Test
	public void testShowDataPoints() {
		testPCA.showDataPoints();
	}
	
	

}
