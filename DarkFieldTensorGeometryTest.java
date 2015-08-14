package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

import edu.stanford.rsl.conrad.geometry.trajectories.Trajectory;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;

public class DarkFieldTensorGeometryTest {

	private int numScatterVectors;
	
	DarkFieldTensorGeometry myGeometry;
	
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		
		Configuration myConfig = new Configuration();
		
		numScatterVectors = 7;
		
		Trajectory geo = myConfig.getGeometry();
		
		geo.setDetectorHeight(100);
		geo.setDetectorWidth(80);
		
		geo.setPixelDimensionX(1f);
		geo.setPixelDimensionY(1f);
		
		geo.setReconDimensionX(50);
		geo.setReconDimensionY(50);
		geo.setReconDimensionZ(50);
		
		geo.setVoxelSpacingX(1d);
		geo.setVoxelSpacingY(1d);
		geo.setVoxelSpacingZ(1d);
		
		geo.setOriginInPixelsX(-25);
		geo.setOriginInPixelsY(-25);
		geo.setOriginInPixelsZ(-25);
		
		geo.setDetectorOffsetU(0);
		geo.setDetectorOffsetV(0);
		
		geo.setAverageAngularIncrement(1);
		
		SimpleVector rotAxis = new SimpleVector(0f,0f,1f);
		geo.setRotationAxis(rotAxis);
		
		myGeometry = new DarkFieldTensorGeometry(myConfig, numScatterVectors);
		
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testCheckEquality() {
		SimpleVector v1 = new SimpleVector(1f,2f,3f);
		SimpleVector v2 = new SimpleVector(2f,2f,3f);
		SimpleVector v3 = new SimpleVector(1f,2f,3f);
		
		double delta = 0.004;
		
		boolean check1 = SimpleOperators.equalElementWise(v1, v2, delta);
		
		assertTrue("checEquality went wrong",check1 == false);
		
		boolean check2 = SimpleOperators.equalElementWise(v1, v3, delta);
		
		assertTrue("checEquality went wrong",check2 == true);
	}
	
	@Test
	public void testCoordinateConversion() {
		
//		double uWorld = 
//		
//		assertTrue("checEquality went wrong",check2 == true);
	}

}
