package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;

public class TensorAbsorptionRecon {

	
	public TensorAbsorptionRecon(){
		
		
	}
	
	public static void main(String[] args){
		
		
		String fileNameConfig1 = "C:\\Users\\schiffers\\workspace\\Configurations\\Config_Full_Resolution_100_cubic.xml";
		// Load configuration wooden case

		Configuration config = Configuration.loadConfiguration(fileNameConfig1);
		System.out.println("Configuration 1 loaded.");
		
		//  Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig2);
		//	System.out.println("Configuration 2 loaded.");

		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		Configuration2.getGeometry().setRotationAxis(rotationAxis2);
				
		// Load ImageJ
		new ImageJ();
		System.out.println("ImageJ started.");

		//Choose dark-field projections. Sinogram is calculated from projections.
		
		String fileNameAMP1 = "E:\\fschiffers\\MeasuredData\\WoodAMP1.tif";
		String fileNameAMP2 = "E:\\fschiffers\\MeasuredData\\WoodAMP2.tif";
		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		ImagePlus imgAMP1 = IJ.openImage(fileNameAMP1);
		DarkField3DSinogram sinoAMP1   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP1);
		
		sinoAMP1.show();
		
		
		// Load dark field image of orientation 2
		//ImagePlus imgAMP2 = IJ.openImage(fileNameAMP2);
		//DarkField3DSinogram sinoAMP2   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP2);		
	
		
		DarkFieldAbsorptionRecon3D reconAMP = new DarkFieldAbsorptionRecon3D(config);
		
		DarkField3DTensorVolume myRecon = reconAMP.reconstructAbsorptionVolume(sinoAMP1);
		
		myRecon.show();

		float th_lower = 0.0005f;
		float th_higher = 0.004f;
		
		DarkField3DTensorVolume myMask = reconAMP.createMaskByBinaryThresholding(th_lower, th_higher);
		
		myMask.show();
		
		
		
		
	}
	
	
}
