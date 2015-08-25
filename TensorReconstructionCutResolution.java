// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;



import java.io.File;

import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume.TensorConstraintType;
import ij.ImageJ;





public class TensorReconstructionCutResolution{

	public static void main (String [] args) throws Exception{

		/*
		 *  INITILIAZE THE DATA FILES
		 */

		// Path to the configuration file (default no second configuration is used)
		String fileNameConfig1 = "E:\\fschiffers\\MeasuredData\\HalfAngleAndCut\\config_CUT.xml";

		// Path to the 2 darkfield images
		
		String fileNameDCI1 = "E:\\fschiffers\\MeasuredData\\HalfAngleAndCut\\WoodDCI1.tif";
		String fileNameDCI2 = "E:\\fschiffers\\MeasuredData\\HalfAngleAndCut\\WoodDCI2.tif";
		
		// Path to the 2 absorption images
		String fileNameAMP1 = "E:\\fschiffers\\MeasuredData\\HalfAngleAndCut\\WoodAMP1.tif";
		String fileNameAMP2 = "E:\\fschiffers\\MeasuredData\\HalfAngleAndCut\\WoodAMP2.tif";
		
		File folder = new File(fileNameDCI1);
		
		/*
		 * INITILIAZATION OF SOME DATA
		 */
		
		float th_lower = 0.0005f;
		float th_higher = 0.004f;
		
		// Number of scatter vectors
		int numScatterVectors = 7;
		//Step size for Gradient decent
		float stepSize = 0.003f;
		// Number of maximal iterations in gradient decent
		int maxIt = 5;
		
		new ImageJ();
		
		// Initialize the pipeline
		DarkFieldReconPipeline myDarkFieldPipeLine = new DarkFieldReconPipeline(fileNameConfig1, fileNameDCI1, null,fileNameConfig1,TensorConstraintType.NO_CONSTRAINT);
		
		// Create the Absorption Mask
		boolean saveAMP = true;
		boolean saveMask = true;
		myDarkFieldPipeLine.reconstructMaskForZeroConstraint(th_lower, th_higher,saveAMP, saveMask,fileNameAMP1);
		
		myDarkFieldPipeLine.getReconMask().show("Mask Image");
		
		
		
		
		System.out.println("Reconstruction mask was successfully created and saved.");
		
		
		boolean writeVtkInEveryStep = true;
		
		// Reconstruct DarkField Volume
		myDarkFieldPipeLine.reconstructDarkFieldVolume(numScatterVectors, maxIt, stepSize, folder,writeVtkInEveryStep);
		
		System.out.println(" DarkField Reconstruction was successfully created and saved.");
		
		myDarkFieldPipeLine.calculateFiberOrientations(true);
		
		System.out.println("Fiber Orientations sucessfully saved.");
		
		
	}

}
