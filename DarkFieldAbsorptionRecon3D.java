package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.ImageToGrid3D;
import edu.stanford.rsl.science.darkfield.parallel.ParallelBackprojectorAMP2D;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;

public class DarkFieldAbsorptionRecon3D  extends  DarkFieldTensorGeometry  {

	DarkField3DTensorVolume reconAMP;
	DarkField3DTensorVolume myMask;
	
	public DarkFieldAbsorptionRecon3D(Configuration config){
		// Call super constructor of TensorGeometry
		super(config,1);
		
	}
	
	
	public DarkField3DTensorVolume reconstructAbsorptionVolume(DarkField3DSinogram sinogram){
		
		
//
//		sinogram.show("Projection Images before Filtering.");
//		sinogram.showSinogram("Sinogram before Filtering.");
		
		RamLakKernel ramLak = new RamLakKernel(maxU_index,deltaU);

		
		
		for(int curTheta = 0; curTheta < maxTheta_index; curTheta++){
		
			Grid2D projImage = sinogram.getSubGrid(curTheta);
			
			for(int curV = 0; curV < maxV_index; curV++){
				// Get Detector row
				Grid1D detectorRow = projImage.getSubGrid(curV);
				// Apply Ram Lak Filter
				ramLak.applyToGrid(detectorRow);
				projImage.setSubGrid(curV, detectorRow);
			}

			sinogram.setSubGrid(curTheta, projImage);
			
		}
		


//		
//		sinogram.show("Projection Images after Filtering");
//		sinogram.showSinogram("Sinogram after Filtering");

		

		 reconAMP = backProjectAbsorption(sinogram);
		
		return reconAMP;
		
	}

	
	public DarkField3DTensorVolume createMask(float th_lower, float th_higher){
		
		if (reconAMP == null){
			return null;
		}
		
		myMask = new DarkField3DTensorVolume(imgSizeX, imgSizeY, imgSizeZ, 1, getSpacing(), getOrigin());
		
		for(int x = 0; x < imgSizeX; x++){
			for(int y = 0; y < imgSizeY; y++){
				for(int z = 0; z < imgSizeZ; z++){
					float val = reconAMP.getAtIndex(x, y, z, 0);
					// If value is in bounds set Mask to 1
					if(val > th_lower && val < th_higher){
						myMask.setValueAtChannelN(x, y, z, 0, 1);
					} else{ // if value is out of bounds set Mask to 0
						myMask.setValueAtChannelN(x, y, z, 0, 0);
					}
				} // END z
				} // END y
				} // END x
		
		return myMask;
	}

	
	public	DarkField3DTensorVolume backProjectAbsorption(DarkField3DSinogram sino3D) {
		
		boolean debug = false;
		
		// Create empty 3DDarkField Volume
		DarkField3DTensorVolume grid = new DarkField3DTensorVolume(imgSizeX,imgSizeY,imgSizeZ,1,getSpacing(),getOrigin());
		
		// Loop over all projection angles
		for (int curTheta = 0; curTheta < maxTheta_index; curTheta++) {

			
			
			// Calculate current projection angle
			double theta = deltaTheta * curTheta;
			
			double cosTheta = Math.cos(theta);
			double sinTheta = Math.sin(theta);

			if(debug){
			// Debug Output: Uncomment if you need to debug the Backprojector
			System.out.println("Cur BackProj: " +curTheta +"/" +maxTheta_index 
					+ " (" +((10000*curTheta/maxTheta_index)/100.0) +"% done.)");
			}
			
			// get detector grid
			Grid2D detectorImageAtTheta = sino3D.getSubGrid(curTheta);		
						
			// Create direction Vector of the detector at given angle Theta
			// Remember: Third coordinate is 0
			SimpleVector dirU = calculateRotatedVector(cosTheta, sinTheta,0).getAbstractVector();
			SimpleVector dirV = calculateRotatedVector(0,0,1f).getAbstractVector();
			
			// Loop through complete volume to do pixel based backprojection
				for (int x = 0; x < imgSizeX; x++) {
					for (int y = 0; y < imgSizeY; y++) {
						for (int z = 0; z < imgSizeZ; z++) {
							
					// compute world coordinate of current pixel
					double[] w = grid.indexToPhysical(x, y,z);
					
					// Create current voxel element
					SimpleVector voxel = new SimpleVector(w[0], w[1],w[2]);
					
					// Calcualte detector coodinates
					SimpleVector orthProj = calcDetectorCoordinates(voxel,dirU,dirV);
					
					// Calculates the subpixel detector coordinate curU
					double curU_index = calcU_index(orthProj.getElement(0));
					// precalculate detector column 						
					double curV_index = calcV_index(orthProj.getElement(1));  
						
					// check detector bounds, continue if out of borders
					if ( 		maxU_index <= curU_index + 1
							||  curU_index < 0
							||  maxV_index < curV_index + 1
							||  curV_index < 0
							){
						continue; // Do nothing if projected point does not lie on detector
					}
					
					// Calculate the interpolated darkField value of the current voxel point
					// of the current projection image
					float ampValue = InterpolationOperators.interpolateLinear(detectorImageAtTheta,curU_index,curV_index);
						
					grid.addAtIndexDarkfield(x, y, z, 0, ampValue);
					
					
					} // END LOOP Z
					} // END LOOP Y
			} // END LOOP X
		} // END LOOP ANGLES
		
		
		
		// TODO Normalization factors
		NumericPointwiseOperators.divideBy(grid, (float) (maxTheta_index / Math.PI));
		return grid;
	}
	
	
	
	
}
