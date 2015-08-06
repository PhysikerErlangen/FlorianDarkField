package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.io.File;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.utils.Configuration;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;

public class ReadFiberDirectionsOutOfTif {

	public static void main(String[] args){

		// Load ImageJ;
		new ImageJ();
		
		
		String pathToDCI_Volume = "E:\\fschiffers\\MeasuredData\\DataTest\\DCI_volume.tif"; 
		
		File fileDCI = new File(pathToDCI_Volume); 
		
		ImagePlus imgDCI1 = IJ.openImage(fileDCI.getPath());
		
		DarkField3DTensorVolume darkFieldVolume = DarkField3DTensorVolume.readFromImagePlus(imgDCI1);
		
		darkFieldVolume.show();
		
		String fileNameConfig1 = "E:\\fschiffers\\Configurations\\Config_Full_Resolution_100_cubic.xml";
		Configuration config1 = Configuration.loadConfiguration(fileNameConfig1);
		DarkFieldScatterCoef scatterCoef = new DarkFieldScatterCoef(config1, 7);
		
		SimpleMatrix scatterMatrix = scatterCoef.getScatterVectorsAsMatrix();
		
		
		DarkFieldReconPipeline.calculateFiberOrientations(true, darkFieldVolume, scatterMatrix, fileDCI);
		
		System.out.println("Everthing was written");
		
	}
	
}
