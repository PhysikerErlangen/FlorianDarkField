package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.ImageToGrid3D;
import edu.stanford.rsl.science.darkfield.parallel.ParallelBackprojectorAMP2D;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;

public class DarkFieldAbsorptionRecon3D  extends  DarkFieldTensorGeometry  {

	public DarkFieldAbsorptionRecon3D(Configuration config, int numScatterVectors){
		// Call super constructor of TensorGeometry
		super(config,numScatterVectors);
		
	}
	
	
	public void reconstructAbsorptionVolume(DarkField3DSinogram absorption){
		
		RamLakKernel ramLak = new RamLakKernel(maxU_index,deltaU);

		
		
		for(int curTheta = 0; curTheta < maxTheta_index; curTheta++){
		
			Grid2D projImage = absorption.getSubGrid(curTheta);
			
			for(int curV = 0; curV < maxU_index; curV++){
				// Get Detector row
				Grid1D detectorRow = projImage.getSubGrid(curV);
				// Apply Ram Lak Filter
				ramLak.applyToGrid(detectorRow);
				projImage.setSubGrid(curV, detectorRow);
			}
		}
		
		
		
	}
	
	// This is a simple Filtered Backprojection Algorithm for the attenuation images
	// The reconstructed result is later used for a zero-constraint which improves the efficieny of Hu's algorithm

	private Grid2D AMPRecon(Grid2D sinogram) {
		
		Grid2D filteredSinogram = new Grid2D(sinogram);
		
		// Create ramLak Kernel
		RamLakKernel ramLak = new RamLakKernel(maxU_index, this.deltaU);
		
		for (int curThetaIndex = 0; curThetaIndex < sinogram.getSize()[1]; ++curThetaIndex) {
		
			ramLak.applyToGrid(filteredSinogram.getSubGrid(curThetaIndex));
	
		}
		
		filteredSinogram.show("The Filtered Sinogram");

		ParallelBackprojectorAMP2D backproj = new ParallelBackprojectorAMP2D(imgSizeX, imgSizeY, (float) spacingX,(float) spacingY, offSetU_index);

		return backproj.backprojectPixelDriven(filteredSinogram);
		
	}


	
	
}
