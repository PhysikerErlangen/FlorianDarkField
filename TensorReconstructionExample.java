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





public class TensorReconstructionExample{

	public static void main (String [] args) throws Exception{

		String fileNameConfig1 = "C:\\Users\\schiffers\\workspace\\Configurations\\Config_Full_Resolution_100_cubic.xml";
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


		//Choose dark-field projections. Sinogram is calculated from projections.
		
		String fileNameDCI1 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\WoodDCI1.tif";
		String fileNameDCI2 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\WoodDCI2.tif";
		
		
		String fileNameAMP1 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\WoodAMP1.tif";
		String fileNameAMP2 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\WoodAMP2.tif";
		
		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		ImagePlus imgDCI1 = IJ.openImage(fileNameDCI1);
		DarkField3DSinogram sinoDCI1   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI1);
		
		// Load dark field image of orientation 2
		ImagePlus imgDCI2 = IJ.openImage(fileNameDCI2);
		DarkField3DSinogram sinoDCI2   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI2);		
	
		// Load absorption image of orientation 1
		ImagePlus imgAMP1 = IJ.openImage(fileNameAMP1);
		DarkField3DSinogram sinoAMP1   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP1);
	
		// Load absorption image of orientation 2
		ImagePlus imgAMP2 = IJ.openImage(fileNameAMP2);
		DarkField3DSinogram sinoAMP2   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP2);		
	
		
boolean showFlag = false;
		
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
	    
		
		imgAMP1 = null;
		imgAMP2 = null;
		imgDCI1 = null;
		imgDCI2 = null;
		
		
		// END loading dark field images 

		// Runtime calculating starts
		long startTime = System.currentTimeMillis();

		
		/*
		 * INITILIAZATION OF SOME DATA
		 */
		
		float th_lower = 0.0005f;
		float th_higher = 0.004f;
		
		System.out.println("Load Parallel Beam Reconstruction Pipeline.");
		// Initialize the Parallel Beam Absorption Reconstruction
		DarkFieldAbsorptionRecon3D parallellBeamRecon = new DarkFieldAbsorptionRecon3D(Configuration1);
		// Reconstruct the Absorption volume later used for zero constraint
		DarkField3DTensorVolume reconAMP = parallellBeamRecon.reconstructAbsorptionVolume(sinoAMP1);
		// Create Mask out of absorption reconstruction
		DarkField3DTensorVolume reconMask = parallellBeamRecon.createMask(th_lower, th_higher);
		if(showFlag)reconMask.show("Mask used for zero constraint in reconstruction.");
		System.out.println("Mask for reconsruction created.");
				
		
		parallellBeamRecon = null;
		
		// Number of scatter vectors
		int numScatterVectors = 3;
		//Stepsize for Gradient decent
		float stepSize = 0.003f;
		// Number of maximal iterations in gradient decent
		int maxIt = 2;
		
		// Initialize the GradientSolver3D

		GradientSolverTensor3D gradientSolver = new GradientSolverTensor3D(Configuration1, Configuration2, sinoDCI1, sinoDCI1, stepSize, maxIt, numScatterVectors,reconMask,reconMask);
		
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
