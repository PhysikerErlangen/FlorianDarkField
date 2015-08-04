// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;



import weka.gui.beans.ConfigurationEvent;

import com.jogamp.opengl.util.awt.ImageUtil;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;


// Used for solving the 3D Gradientsolver in the tensor framework
import edu.stanford.rsl.science.darkfield.FlorianDarkField.GradientSolverTensor3D;

// Contains the reconstructed sample
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ImageToSinogram3D;
import edu.stanford.rsl.science.darkfield.parallel.ParallelBackprojectorAMP2D;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;





public class TensorReconstructionWoodLowResolution{

	public static void main (String [] args) throws Exception{

		String fileNameConfig1 = "C:\\Users\\schiffers\\workspace\\Configurations\\HalfAngleCropped_Cut_60_38_100_Date_07-27-15.xml";
		
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
		
		String fileNameDCI1 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\Reduced_Data_Cut\\HalfAngleCropped_WoodDCI1.tif";
		String fileNameDCI2 = "C:\\Users\\schiffers\\workspace\\MeasuredData\\Reduced_Data_Cut\\HalfAngleCropped_WoodDCI2.tif";
		
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
		int numScatterVectors = 7;
		//Stepsize for Gradient decent
		float stepSize = 0.01f;
		// Number of maximal iterations in gradient decent
		int maxIt = 1;
		
		// Initialize the GradientSolver3D
		GradientSolverTensor3D gradientSolver = new GradientSolverTensor3D(Configuration1, Configuration2, sinoDCI1, sinoDCI2, stepSize, maxIt, numScatterVectors);
		
		DarkField3DTensorVolume reconImage = gradientSolver.Gradient3D();

		
		ImagePlus test = DarkField3DTensorVolume.wrapDarkFieldGrid3DTensorToImagePlus(reconImage, "test");
		
		test.show();
		
		IJ.save(test,"C:\\Users\\schiffers\\workspace\\MeasuredData\\testObjectRecon.tif");
		
       //reconImage.show();
		
		
		/*
		 * END LOAD CONIFUGRATION PARAMETERS
		 */
		
		long endTime = System.currentTimeMillis();
		System.out.println("Whole tensor reconstruction was done in " +(endTime-startTime) + "ms.");

	}

	

	
}
