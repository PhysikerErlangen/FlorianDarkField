/**
 * 
 */
package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;
import ij.ImageJ;

import java.io.File;
import java.util.Scanner;

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

	DarkFieldTensorPhantom phantomObject;

	DarkFieldReconPipeline myDarkFieldPipeLine;
	
	int numScatterVectors;
	
	File folder;
	
	/**
	 * @throws java.lang.Exception
	 */
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		  Scanner sc = new Scanner(System.in);
		    System.out.print("Gib was ein: ");
		    String eingabe = sc.next();
		    System.out.println("Du hast " + eingabe + " eingegeben.");
		    sc.close();
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		
	String fileNameConfig1 = "E:\\fschiffers\\MeasuredData\\Phantom2\\PhantomHalfLarge_unsymetricSmall.xml";
		
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
		phantomObject = new DarkFieldTensorPhantom(Configuration1,numScatterVectors);

		// display the phantom
		phantomObject.getPhantom().show("Phantom Volume");
		
		phantomObject.calculateDarkFieldProjection(Configuration1, Configuration2);
 		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		DarkField3DSinogram sinoDCI1   = phantomObject.getDarkFieldSinogram(0);		
		DarkField3DSinogram sinoDCI2   = phantomObject.getDarkFieldSinogram(1);
		
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
		float stepSize = 0.01f;
		// Number of maximal iterations in gradient decent
		int maxIt =20;
		
		// Initialize the pipeline
		myDarkFieldPipeLine = new DarkFieldReconPipeline(Configuration1,Configuration2,fileNameConfig1);
		
		// Reconstruct DarkField Volume
		
		folder = new File(fileNameConfig1);
		
		myDarkFieldPipeLine.reconstructDarkFieldVolume(numScatterVectors,maxIt,stepSize,folder,sinoDCI1,sinoDCI2);
		
		System.out.println(" DarkField Reconstruction was successfully created and saved.");
		
		File myParentFile = new File(fileNameConfig1);
		
		myDarkFieldPipeLine.calculateFiberOrientations(myParentFile);
		
		System.out.println("Fiber Orientations sucessfully saved.");
		
	}

	// This methods tests, if the creation of the Phantom works correctly
	@Test
	public void testPhantomReconDifference() {
		
		DarkField3DTensorVolume recon =  myDarkFieldPipeLine.getReconDarkField();
		
		DarkField3DTensorVolume  diffVolume  = DarkField3DTensorVolume.sub(recon, phantomObject.getPhantom());
		
		diffVolume.show("Difference between diffVolume and Phantom");
		
		SimpleMatrix myScatterMatrix = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		DarkFieldReconPipeline.calculateFiberOrientations(phantomObject.getPhantom(), myScatterMatrix, folder,"fiberDirectionsPhantom");
		
	}

}
