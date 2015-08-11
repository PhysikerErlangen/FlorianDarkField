/**
 * 
 */
package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;
import ij.ImageJ;

import java.io.File;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;

/**
 * @author schiffers
 *
 */
public class DarkFieldTensorPhantomTest {

	DarkFieldTensorPhantom phantom;

	DarkFieldReconPipeline myDarkFieldPipeLine;
	
	int numScatterVectors;
	
	File folder;
	
	/**
	 * @throws java.lang.Exception
	 */
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		
	String fileNameConfig1 = "E:\\fschiffers\\MeasuredData\\Phantom2\\PhantomHalfLarge_unsymetric.xml";
		
		// Load configuration wooden case

		Configuration Configuration1 = Configuration.loadConfiguration(fileNameConfig1);
		System.out.println("Configuration 1 loaded.");

		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		Configuration2.getGeometry().setRotationAxis(rotationAxis2);
		
		
		// Load ImageJ
		new ImageJ();
		
		// Number of scatter vectors
		numScatterVectors = 7;
		
		// Create Dark Field Phantom
		phantom = new DarkFieldTensorPhantom(Configuration1,numScatterVectors);

		// display the phantom
		phantom.phantom.show("Phantom Volume");
		
		phantom.calculateDarkFieldProjection(Configuration1, Configuration2);
 		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		DarkField3DSinogram sinoDCI1   = phantom.getDarkFieldSinogram(0);		
		
		boolean showFlag = true;
		
		if(showFlag){
		
		// Show DarkField Projection Images
		sinoDCI1.show("Dark field image 1");
		// Show Stack of DarkField slices
		sinoDCI1.showSinogram("Sinogram of dark field image 1");
		// Show DarkField Projection Images
		}

		
		/*
		 * INITILIAZATION OF SOME DATA
		 */
		
		// Number of scatter vectors
		//Step size for Gradient decent
		float stepSize = 0.007f;
		// Number of maximal iterations in gradient decent
		int maxIt = 1;
		
		// Initialize the pipeline
		myDarkFieldPipeLine = new DarkFieldReconPipeline(Configuration1,Configuration2);
		
		// Reconstruct DarkField Volume
		
		folder = new File(fileNameConfig1);
		
		myDarkFieldPipeLine.reconstructDarkFieldVolume(numScatterVectors,maxIt,stepSize,folder,sinoDCI1,null);
		
		System.out.println(" DarkField Reconstruction was successfully created and saved.");
		
		File myParentFile = new File(fileNameConfig1);
		
		myDarkFieldPipeLine.calculateFiberOrientations(myParentFile);
		
		System.out.println("Fiber Orientations sucessfully saved.");
		
	}


	// This methods tests, if the creation of the Phantom works correctly
	@Test
	public void testPhantomFiberExtraction() {
		
		SimpleMatrix myScatterMatrix = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		DarkFieldReconPipeline.calculateFiberOrientations(phantom.phantom, myScatterMatrix, folder,"fiberDirectionsPhantom");
		
		
		fail("Not yet implemented");
	}

}
