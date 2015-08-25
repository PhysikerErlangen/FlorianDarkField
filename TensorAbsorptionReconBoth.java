package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;

public class TensorAbsorptionReconBoth {

	
	public TensorAbsorptionReconBoth(){
		
		
	}
	
	public static void main(String[] args){
		
		
		String fileNameConfig1 = "E:\\fschiffers\\Configurations\\Config_Full_Resolution_100_cubic.xml";
		// Load configuration wooden case

		Configuration Config1 = Configuration.loadConfiguration(fileNameConfig1);
		System.out.println("Configuration 1 loaded.");
		
		//  Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig2);
		//	System.out.println("Configuration 2 loaded.");

		Configuration Config2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		Config2.getGeometry().setRotationAxis(rotationAxis2);
				
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
		
		sinoAMP1.show("AMP Signal 1");
		
		
		// Load dark field image of orientation 2
		ImagePlus imgAMP2 = IJ.openImage(fileNameAMP2);
		DarkField3DSinogram sinoAMP2   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP2);		
		
		sinoAMP2.show("AMP Signal 2");

		/*
		 * Reconstruction
		 */
		
		DarkFieldAbsorptionRecon3D reconAMP1 = new DarkFieldAbsorptionRecon3D(Config1);
		DarkField3DTensorVolume myRecon1 = reconAMP1.reconstructAbsorptionVolume(sinoAMP1);
		myRecon1.show("Recon AMP 1");
		
		DarkFieldAbsorptionRecon3D reconAMP2 = new DarkFieldAbsorptionRecon3D(Config2);
		DarkField3DTensorVolume myRecon2 = reconAMP2.reconstructAbsorptionVolume(sinoAMP2);
		myRecon2.show("Recon AMP 2");
		
		DarkField3DTensorVolume addRecon = DarkField3DTensorVolume.add(myRecon1,myRecon2);
		addRecon.show("Recon addition both");

		

//		float th_lower = 0.0005f;
//		float th_higher = 0.004f;
//		
//		DarkField3DTensorVolume myMask = reconAMP.createMaskByBinaryThresholding(th_lower, th_higher);
//		
//		myMask.show();
		
		
		
		
	}
	
	
}
