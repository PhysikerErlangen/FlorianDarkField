// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;



import com.jogamp.opengl.util.awt.ImageUtil;

import edu.stanford.rsl.conrad.utils.Configuration;


// Used for solving the 3D Gradientsolver in the tensor framework
import edu.stanford.rsl.science.darkfield.FlorianDarkField.GradientSolverTensor3D;

// Contains the reconstructed sample
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ImageToSinogram3D;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;





public class TensorReconstructionExample{

	public static void main (String [] args) throws Exception{

		String fileNameConfig1 = "C:\\Users\\schiffers\\workspace\\Configurations\\ConfigurationFlorian-001-axis.xml";
		String fileNameConfig2 = "C:\\Users\\schiffers\\workspace\\Configurations\\ConfigurationFlorian-010-axis.xml";
		
		// Load configuration wooden case

		Configuration Configuration1 = Configuration.loadConfiguration(fileNameConfig1);
		System.out.println("Configuration 1 loaded.");
		
		
		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig2);
		System.out.println("Configuration 2 loaded.");
		// Load ImageJ
		new ImageJ();
		System.out.println("ImageJ started.");


		//Choose dark-field projections. Sinogram is calculated from projections.
		
		String fileNameDCI1 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\WoodDCI1.tif";
		String fileNameDCI2 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\WoodDCI2.tif";
		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		ImagePlus imgDCI1 = IJ.openImage(fileNameDCI1);
		DarkField3DSinogram sinoDCI1   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI1);
		
		// Load dark field image of orientation 2
		ImagePlus imgDCI2 = IJ.openImage(fileNameDCI2);
		DarkField3DSinogram sinoDCI2   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI2);		
	

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
	    
		// END loading dark field images 

		// Runtime calculating starts
		long startTime = System.currentTimeMillis();

		
		/*
		 * INITILIAZATION OF SOME DATA
		 */
		
		// Number of scatter vectors
		int numScatterVectors = 3;
		//Stepsize for Gradient decent
		float stepSize = 0.003f;
		// Number of maximal iterations in gradient decent
		int maxIt = 2;
		
		// Initialize the GradientSolver3D
		GradientSolverTensor3D gradientSolver = new GradientSolverTensor3D(Configuration1, Configuration2, sinoDCI1, sinoDCI1, stepSize, maxIt, numScatterVectors);
		
		DarkField3DTensorVolume reconImage = gradientSolver.Gradient3D();

		
		ImagePlus test = edu.stanford.rsl.conrad.utils.ImageUtil.wrapGrid4D(reconImage.getMultichannelData(), "test");
		
		IJ.save(test,"C:\\Users\\schiffers\\workspace\\MeasuredData\\testObjectRecon.tif");
		
       reconImage.show();
		
		
		/*
		 * END LOAD CONIFUGRATION PARAMETERS
		 */
		
		long endTime = System.currentTimeMillis();
		System.out.println("Whole tensor reconstruction was done in " +(endTime-startTime) + "ms.");

	}

}
