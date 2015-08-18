package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.ImageToGrid3D;

public class TensorAbsorptionRecon2Traj {

	
	public TensorAbsorptionRecon2Traj(){
		
		
	}
	
	public static void main(String[] args){
		
		
		String fileNameConfig = "E:\\fschiffers\\Configurations\\Config_Full_Resolution_010_cubic.xml";
		// Load configuration wooden case

		Configuration config = Configuration.loadConfiguration(fileNameConfig);
		System.out.println("Configuration loaded.");
		
		//  Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig2);
		//	System.out.println("Configuration 2 loaded.");

		
		SimpleVector myRotationAxis = config.getGeometry().getRotationAxis();
		System.out.println("Rotation axis is: " + myRotationAxis);
		
		double offSetU = config.getGeometry().getDetectorOffsetU();
		System.out.println("OffsetU is: " + offSetU);
		double offSetV = config.getGeometry().getDetectorOffsetV();
		System.out.println("OffsetV is: " + offSetV);
				
		// Load ImageJ
		new ImageJ();
		System.out.println("ImageJ started.");

		//Choose dark-field projections. Sinogram is calculated from projections.
		
		String fleNameAMP1 = "E:\\fschiffers\\MeasuredData\\WoodAMP2.tif";
		// String fileNameAMP2 = "E:\\fschiffers\\MeasuredData\\WoodAMP2.tif";
		
		/* 
		 * Load dark field images
		 */
		
		// Load dark field image of orientation 1
		ImagePlus imgAMP = IJ.openImage(fleNameAMP1);
		DarkField3DSinogram sinoAMP   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP);
		
		sinoAMP.show();
		
		
		// Load dark field image of orientation 2
		//ImagePlus imgAMP2 = IJ.openImage(fileNameAMP2);
		//DarkField3DSinogram sinoAMP2   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP2);		
	
		
		DarkFieldAbsorptionRecon3D reconAMP = new DarkFieldAbsorptionRecon3D(config);
		
		DarkField3DTensorVolume myRecon = reconAMP.reconstructAbsorptionVolume(sinoAMP);
		
		myRecon.show();

		float th_lower = 0.0005f;
		float th_higher = 0.004f;
		
		DarkField3DTensorVolume myMask = reconAMP.createMaskByBinaryThresholding(th_lower, th_higher);
		
		myMask.show();
		
		
		
		
	}
	
	
}
