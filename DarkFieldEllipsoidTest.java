package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldEllipsoidTest {

	
	DarkFieldEllipsoid myEllipsoid;
	
	int numScatterVectors;
	
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

		numScatterVectors = 7;
		
		// Get scatterDirections
		SimpleMatrix scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		// Preset CoefInit to 1. This is just for a "starting value" to calculate projections
		SimpleVector scatterCoef = new SimpleVector(numScatterVectors);
		for(int channel = 0; channel < numScatterVectors; channel ++){
			scatterCoef.setElementValue(channel, 1);
		}
	
		// Preset eigenValues
		double lambda1 = 3;
		double lambda2 = 3;
		double lambda3 = 3;
		SimpleVector eigenValues  = new SimpleVector(lambda1,lambda2,lambda3);
		
		// Preset EigenVectors
		SimpleVector v1 = new SimpleVector(0,1,0);
		SimpleVector v2 = new SimpleVector(0,0,1);
		SimpleVector v3 = new SimpleVector(1,0,0);
		SimpleMatrix eigenVectors = new SimpleMatrix(3,3);
		eigenVectors.setColValue(0, v1);
		eigenVectors.setColValue(1, v2);
		eigenVectors.setColValue(2, v3);
	
		myEllipsoid = new DarkFieldEllipsoid(scatterDirections, scatterCoef, eigenValues, eigenVectors);
		
	}
	

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void TesttransformPointIntoEigenVectorSystem() {
		
		SimpleVector eX = new SimpleVector(1,0,0);
		SimpleVector eY = new SimpleVector(0,1,0);
		SimpleVector eZ = new SimpleVector(0,0,1);
		
		
		SimpleVector point1 = new SimpleVector(1,1,0).normalizedL2();
		SimpleVector point2 = new SimpleVector(-1,1,0).normalizedL2();
		SimpleVector point3 = new SimpleVector(0,0,1).normalizedL2();
		
		SimpleVector vec1 = new SimpleVector(1,1,0).normalizedL2();
		SimpleVector vec2 = new SimpleVector(-1,1,0).normalizedL2();
		SimpleVector vec3 = new SimpleVector(0,0,1).normalizedL2();
		
		SimpleMatrix eigenVectors= new SimpleMatrix(3,3);
		eigenVectors.setColValue(0, vec1);
		eigenVectors.setColValue(1, vec2);
		eigenVectors.setColValue(2, vec3);
		
		SimpleVector transX  = DarkFieldEllipsoid.transformPointIntoEigenVectorSystem(point1, eigenVectors);
		SimpleVector transY  = DarkFieldEllipsoid.transformPointIntoEigenVectorSystem(point2, eigenVectors);
		SimpleVector transZ  = DarkFieldEllipsoid.transformPointIntoEigenVectorSystem(point3, eigenVectors);

		double delta = 1E-5;
		
		assertTrue("Not equal", SimpleOperators.equalElementWise(transX,eX, delta));
		assertTrue("Not equal", SimpleOperators.equalElementWise(transY,eY, delta));
		assertTrue("Not equal", SimpleOperators.equalElementWise(transZ,eZ, delta));
		
	}

	
	@Test
	public void testCalculateProjectedCoefficients() {
		
		SimpleVector projScatterCoef = myEllipsoid.calculateSquaredProjectedCoefficients();
		
		System.out.println("Projected Scatter Coefficients:");
		System.out.println(projScatterCoef);
		
		
		
		
	}

}
