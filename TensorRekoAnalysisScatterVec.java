// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;



import java.io.File;
import java.util.ArrayList;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;


import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkFieldTensorPhantom.PhantomType;
// Contains the reconstructed sample
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume.TensorConstraintType;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkFieldTensorGeometry.TrajectoryType;
import ij.ImageJ;



public class TensorRekoAnalysisScatterVec{

	String fileNameConfig1;
	
	DarkField3DSinogram sinoDCI1;
	DarkField3DSinogram sinoDCI2;
	
	DarkField3DTensorVolume mask;
	
	DarkFieldTensorClass phantomDirection;
	
	public static void main (String [] args) throws Exception{

		String directory = "E:\\fschiffers\\MeasuredData\\PhantomAnalysis\\";
		String fileName = "PhantomHalfLarge_unsymetricSmall.xml";

		TensorRekoAnalysisScatterVec myAnalysis = new TensorRekoAnalysisScatterVec(directory,fileName);
		
	}
	
	
	public TensorRekoAnalysisScatterVec(String directory,String fileName){

		// Load ImageJ
		new ImageJ();
		
		this.fileNameConfig1 = directory + fileName;
		
		int numScatterDirectionsPhantom  = 100;
		PhantomType phantomType = PhantomType.WOODEN_BLOCK_PHANTOM;
		// Create Phantom data projections
		createProjections(numScatterDirectionsPhantom,phantomType);
		
		DarkFieldReconPipeline myDarkFieldPipeLine;
		
		
		File folder;
		

		// Load configuration wooden case

		Configuration Configuration1 = Configuration.loadConfiguration(fileNameConfig1);
		System.out.println("Configuration 1 loaded.");

		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		Configuration2.getGeometry().setRotationAxis(rotationAxis2);
		System.out.println("Configuration 2 loaded.");
		
		
		// Initialize the pipeline
		myDarkFieldPipeLine = new DarkFieldReconPipeline(Configuration1,Configuration2,fileNameConfig1,TensorConstraintType.NO_CONSTRAINT);
		myDarkFieldPipeLine.setReconMask(mask);
		
		
		/*
		 * INITILIAZATION OF SOME DATA
		 */
		
		// Number of scatter vectors
		//Step size for Gradient decent
		float stepSize = 0.02f;
		// Number of maximal iterations in gradient decent
		int maxIt = 30;
		

		boolean writeVtkInEveryStep = false;
		

		int startNum = 8;
		int endNum  = 30;
		
		SimpleMatrix angularErrorMatrix = new SimpleMatrix(endNum - startNum,2);
		SimpleMatrix sinoDiffErrorMatrix = new SimpleMatrix(endNum - startNum,maxIt);
		
		
		
		for(int numScatterVectors = startNum; numScatterVectors < endNum; numScatterVectors++){
		
		/*
		 * Create new folder
		 */
			folder = new File(directory + "scatterDir" + numScatterVectors);
			folder.mkdir();
			folder = new File(directory + "scatterDir" + numScatterVectors + "\\test.tif");
		
		
		SimpleMatrix errorDiffSino = myDarkFieldPipeLine.reconstructDarkFieldVolume(numScatterVectors,maxIt,stepSize,folder,sinoDCI1,sinoDCI2,writeVtkInEveryStep);
		
		
		
		System.out.println(" DarkField Reconstruction was successfully created and saved.");
		
		File myParentFile = new File(fileNameConfig1);
		
		DarkFieldTensorClass rekoFiberDirection = myDarkFieldPipeLine.calculateFiberOrientations(myParentFile);
		
		double angularError = DarkFieldErrorMeasures.errorAngularDistance(phantomDirection,rekoFiberDirection,mask);
		
		angularErrorMatrix.setRowValue(numScatterVectors-startNum,new SimpleVector(numScatterVectors,angularError));
		sinoDiffErrorMatrix.setRowValue(numScatterVectors-startNum, errorDiffSino.getCol(1));
		
		System.out.println("Fiber Orientations sucessfully saved.");
		
		//DarkField3DTensorVolume recon =  myDarkFieldPipeLine.getReconDarkField();
		
		}
		
		DarkFieldErrorMeasures.writeErrorToTxt(new File(fileNameConfig1) , "angularErrorDifferentScatterVectors.txt", angularErrorMatrix);
		DarkFieldErrorMeasures.writeErrorToTxt(new File(fileNameConfig1) , "SinoDifferennce.txt", sinoDiffErrorMatrix);
		
	}
	
	
	public void createProjections(int numScatterVectors, PhantomType phantomType){
		
		DarkFieldTensorPhantom phantomObject;
		
		
		

		// Load configuration wooden case

		Configuration Configuration1 = Configuration.loadConfiguration(fileNameConfig1);
		System.out.println("Configuration 1 loaded.");

		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		Configuration2.getGeometry().setRotationAxis(rotationAxis2);
	
		
		
		
		// Create Dark Field Phantom
		phantomObject = new DarkFieldTensorPhantom(Configuration1,numScatterVectors,phantomType);

		// display the phantom
		phantomObject.getPhantom().show("Phantom Volume");
		
		File folder = new File(fileNameConfig1);
		SimpleMatrix myScatterMatrix = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		phantomDirection = DarkFieldReconPipeline.calculateFiberOrientations(phantomObject.getPhantom(), myScatterMatrix, folder,"fiberDirectionsPhantom");
		
		System.out.println("Start calculating Phantom DarkField Projections");
		phantomObject.calculateDarkFieldProjection(Configuration1, Configuration2);
 		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		sinoDCI1   = phantomObject.getDarkFieldSinogram(TrajectoryType.HORIZONTAL);
		sinoDCI1.writeDarkFieldSinogramToImage(folder, "DarkFieldSinoGram1.tif","DarkField Sinogram 1");
		sinoDCI2   = phantomObject.getDarkFieldSinogram(TrajectoryType.VERTICAL);
		sinoDCI2.writeDarkFieldSinogramToImage(folder, "DarkFieldSinoGram2.tif","DarkField Sinogram 2");
		
		boolean showFlag = true;
		
		if(showFlag){
		
		// Show DarkField Projection Images
		sinoDCI1.show("Dark field image 1");
		// Show Stack of DarkField slices
		sinoDCI1.showSinogram("Sinogram of dark field image 1");
		// Show DarkField Projection Images
		
		// Show DarkField Projection Images
		sinoDCI2.show("Dark field image 2");
		// Show Stack of DarkField slices
		sinoDCI2.showSinogram("Sinogram of dark field image 2");
		// Show DarkField Projection Images

		
		
		}

		mask = phantomObject.getPhantomMask();
		// mask.show("Mask volume");
		
		
		
		
	}

}
