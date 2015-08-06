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





public class TensorReconstructionPhantom{

	public static void main (String [] args) throws Exception{
 
		String fileNameConfig1 = "E:\\fschiffers\\Configurations\\PhantomHalfLarge_unsymetric.xml";
		
		// Load configuration wooden case

		Configuration Configuration1 = Configuration.loadConfiguration(fileNameConfig1);
		System.out.println("Configuration 1 loaded.");
		
//		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig2);
//		System.out.println("Configuration 2 loaded.");

		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		Configuration2.getGeometry().setRotationAxis(rotationAxis2);
		
		
		
		// Load ImageJ
		new ImageJ();
		System.out.println("ImageJ started.");
		
		
		// Number of scatter vectors
		int numScatterVectors = 3;


		
		// Create Dark Field Phantom
		DarkFieldTensorPhantom phantom = new DarkFieldTensorPhantom(Configuration1,numScatterVectors);
		

		// display the phantom
		phantom.phantom.show("Phantom Volume");
		
		phantom.calculateDarkFieldProjection(Configuration1, Configuration2);
 		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		DarkField3DSinogram sinoDCI1   = phantom.getDarkFieldSinogram(0);		
		
		
		// Load dark field image of orientation 2
		DarkField3DSinogram sinoDCI2   = phantom.getDarkFieldSinogram(1);		


		boolean showFlag = true;
		
		if(showFlag){
		
		// Show DarkField Projection Images
		sinoDCI1.show("Dark field image 1");
		// Show Stack of DarkField slices
		sinoDCI1.showSinogram("Sinogram of dark field image 1");
		// Show DarkField Projection Images
		sinoDCI2.show("Dark field image 2");
		// Show Stack of DarkField slices
	    sinoDCI2.showSinogram("Sinogram of dark field image 2");

		}

		
		/*
		 * INITILIAZATION OF SOME DATA
		 */
		
		
		//Stepsize for Gradient decent
		float stepSize = 0.008f;
		// Number of maximal iterations in gradient decent
		int maxIt = 10;
		
		// Initialize the GradientSolver3D
		GradientSolverTensor3D gradientSolver = new GradientSolverTensor3D(Configuration1, Configuration2, sinoDCI1, sinoDCI1, stepSize, maxIt, numScatterVectors);
		
		DarkField3DTensorVolume reconImage = gradientSolver.Gradient3D();
		
		reconImage.show();
		
		System.out.println("END ALGORITHM");
		//ImagePlus test = edu.stanford.rsl.conrad.utils.ImageUtil.wrapGrid4D(reconImage.getMultichannelData(), "test");
		
//		IJ.save(test,"C:\\Users\\schiffers\\workspace\\MeasuredData\\testObjectRecon.tif");
		
       
//		
//		
//		/*
//		 * END LOAD CONIFUGRATION PARAMETERS
//		 */
//		
//		long endTime = System.currentTimeMillis();
//		System.out.println("Whole tensor reconstruction was done in " +(endTime-startTime) + "ms.");

	}

}
