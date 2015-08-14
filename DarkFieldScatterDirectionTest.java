package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;

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
		
		int numScatterVectors = 7;
		
		SimpleMatrix scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		SimpleVector scatterCoef = new SimpleVector(numScatterVectors);
		for(int channel = 0; channel < numScatterVectors; channel++){
			scatterCoef.setElementValue(channel, 1);
		}
		
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterDirections, scatterCoef);
		
		myPCA.showDataPoints();
		

	}

	@Test
public void testAnisotropic() {
		
		int numScatterVectors = 7;
		
		SimpleMatrix scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		SimpleVector scatterCoef = new SimpleVector(1f,2f,3f,4f,5f,6f,7f);
		
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterDirections, scatterCoef);
		
		myPCA.showDataPoints();
		
	}
	
}
