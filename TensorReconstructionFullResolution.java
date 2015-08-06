// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;



import com.jogamp.opengl.util.awt.ImageUtil;

import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;


// Used for solving the 3D Gradientsolver in the tensor framework
import edu.stanford.rsl.science.darkfield.FlorianDarkField.GradientSolverTensor3D;

// Contains the reconstructed sample
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ImageToSinogram3D;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;





public class TensorReconstructionFullResolution{

	public static void main (String [] args) throws Exception{

		/*
		 *  INITILIAZE THE DATA FILES
		 */

		// Path to the configuration file (default no second configuration is used)
		String fileNameConfig1 = "E:\\fschiffers\\Configurations\\Config_Full_Resolution_100_cubic.xml";

		// Path to the 2 darkfield images
		
		String fileNameDCI1 = "E:\\fschiffers\\MeasuredData\\WoodDCI1.tif";
		String fileNameDCI2 = "E:\\fschiffers\\MeasuredData\\WoodDCI2.tif";
		
		// Path to the 2 absorption images
		String fileNameAMP1 = "E:\\fschiffers\\MeasuredData\\WoodAMP1.tif";
		String fileNameAMP2 = "E:\\fschiffers\\MeasuredData\\WoodAMP2.tif";
		
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
		
		// Initialize the pipeline
		DarkFieldReconPipeline myDarkFieldPipeLine = new DarkFieldReconPipeline(fileNameConfig1, fileNameDCI1, fileNameDCI2, fileNameAMP1, fileNameAMP2);
		
		// Create the Absorption Mask
		boolean saveAMP = true;
		boolean saveMask = true;
		myDarkFieldPipeLine.reconstructMaskForZeroConstraint(th_lower, th_higher,saveAMP, saveMask);
		
		System.out.println("Reconstruction mask was successfully created and saved.");
		
		// Reconstruct DarkField Volume
		boolean saveDarkField = true;
		myDarkFieldPipeLine.reconstructDarkFieldVolume(numScatterVectors, maxIt, stepSize, saveDarkField);
		
		System.out.println(" DarkField Reconstruction was successfully created and saved.");
		
		myDarkFieldPipeLine.calculateFiberOrientations(true);
		
		System.out.println("Fiber Orientations sucessfully saved.");
		
		
	}

}
