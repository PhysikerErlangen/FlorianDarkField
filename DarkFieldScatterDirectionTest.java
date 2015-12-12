package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.util.Scanner;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldScatterDirectionTest {

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		
		  Scanner sc = new Scanner(System.in);
		    System.out.print("Insert anything to close test.");
		    String input = sc.next();
		    System.out.println("You've input: " +input);
		    sc.close();
		
	}

	
	
	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testIsotropicCoefficients() {
		
		for(int k = 25; k < 26; k++){
		
		SimpleMatrix scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(k);
		
		SimpleVector scatterCoef = new SimpleVector(k);
		for(int channel = 0; channel < k; channel++){
			scatterCoef.setElementValue(channel, 1);
			
		}
		
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterDirections, scatterCoef);
		
		myPCA.showDataPoints();
		
		}
		

	}


//	@Test
//public void testAnisotropic7() {
//		
//		int numScatterVectors = 7;
//		
//		SimpleMatrix scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
//		
//		SimpleVector scatterCoef = new SimpleVector(1f,2f,3f,4f,5f,6f,7f);
//		
//		DarkFieldPCA myPCA = new DarkFieldPCA(scatterDirections, scatterCoef);
//		
//		myPCA.showDataPoints();
//		
//	}
//	
//	@Test
//public void testAnisotropic13() {
//		
//		int numScatterVectors = 13;
//		
//		SimpleMatrix scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
//		
//		SimpleVector scatterCoef = new SimpleVector(1f,2f,3f,4f,5f,6f,7f,8f,9f,10f,11f,12f,13f);
//		
//		DarkFieldPCA myPCA = new DarkFieldPCA(scatterDirections, scatterCoef);
//		
//		myPCA.showDataPoints();
//		
//	}
	
}
