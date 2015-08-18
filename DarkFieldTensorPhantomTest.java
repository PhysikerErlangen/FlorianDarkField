package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;
import ij.ImageJ;

import java.util.Scanner;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

import edu.stanford.rsl.conrad.geometry.trajectories.Trajectory;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;

public class DarkFieldTensorPhantomTest {

	DarkFieldTensorPhantom myPhantom;
	
	DarkFieldTensorPhantom myMask;
	
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
		
		Configuration myConfig = new Configuration();
		
		int numScatterVectors = 7;
		
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
		
		
		myPhantom = new DarkFieldTensorPhantom(myConfig, numScatterVectors);
		
		
		
		myPhantom.calcPhantom();
		
		new ImageJ();
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testShowPhantom() {
		myPhantom.getPhantom().show("Phantom Volume");
	}
	
	@Test
	public void testFiberDirection1() {

		SimpleVector fiberDir1 = new SimpleVector(1f,0f,0f);
		System.out.println("Fiber Direction:");
		System.out.println(fiberDir1);
		SimpleVector scatterCoef = DarkFieldTensorPhantom.calculateEllipsoidFromFiberOrientation(fiberDir1,myPhantom.getScatterDirections());
		System.out.println("Projected scatter coefficients:");
		System.out.println(scatterCoef);
		
		SimpleMatrix myScatterTensorPoints = DarkFieldPCA.calcSetOfScatterPoints(myPhantom.getScatterDirections(),scatterCoef);
		
		DarkFieldPointCloudViewer myPointViewer = new DarkFieldPointCloudViewer(myScatterTensorPoints);
		myPointViewer.showPoints();
		
		
		}

	@Test
	public void testFiberDirection2() {

		SimpleVector fiberDir2 = new SimpleVector(1f,0f,0f);
		System.out.println("Fiber Direction:");
		System.out.println(fiberDir2);
		SimpleVector scatterCoef = DarkFieldTensorPhantom.calculateEllipsoidFromFiberOrientation(fiberDir2,myPhantom.getScatterDirections());
		System.out.println("Projected scatter coefficients:");
		System.out.println(scatterCoef);
		
		SimpleMatrix myScatterTensorPoints = DarkFieldPCA.calcSetOfScatterPoints(myPhantom.getScatterDirections(),scatterCoef);
		
		DarkFieldPointCloudViewer myPointViewer = new DarkFieldPointCloudViewer(myScatterTensorPoints);
		myPointViewer.showPoints();
		
		
		}

	
}
