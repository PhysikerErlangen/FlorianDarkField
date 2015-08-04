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

	private DarkField3DTensorVolume reconAMP;
	
	public DarkField3DTensorVolume getReconAMP() {
		return reconAMP;
	}


	public void setReconAMP(DarkField3DTensorVolume reconAMP) {
		this.reconAMP = reconAMP;
	}


	private DarkField3DTensorVolume myMask;
	
	public DarkField3DTensorVolume getMyMask() {
		return myMask;
	}


	public void setMyMask(DarkField3DTensorVolume myMask) {
		this.myMask = myMask;
	}


	public DarkFieldAbsorptionRecon3D(Configuration config){
		// Call super constructor of TensorGeometry
		super(config,1);
		
	}
	
	
	/**
	 * Reconstructs the 3D Absorption Volume by a given Sinogram
	 * @param sinogram - DarkField3DSinogram class is used for convenience
	 * even though the sinogram is a standard absorption sinogram. But can be used for both
	 * @return - Reconstructed Absorption volume
	 */
	public DarkField3DTensorVolume reconstructAbsorptionVolume(DarkField3DSinogram sinogram){
		// Create RamLak Kernel
		RamLakKernel ramLak = new RamLakKernel(maxU_index,deltaU);
		// Loop through all projections
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
		
		// Backprojected filtered sinogram
		 reconAMP = backProjectAbsorption(sinogram);
		
		return reconAMP;
		
	}

	
	/**
	 * Creates a mask by binary thresholding of the reconstructed absorption volume
	 * @param th_lower - lower threshold
	 * @param th_higher - higher thresholds
	 * @return - Volume Mask with 1's and 0'2.
	 * Even though it's an Mask and could be represented as an Grid3D we use 
	 * DarkField3DTensorVolume out of convenience
	 */
	public DarkField3DTensorVolume createMaskByBinaryThresholding(float th_lower, float th_higher){
		
		if (reconAMP == null){
			return null;
		}
		
		// Allocate an Absorption volume
		// Caution: NumScatterVectors set to 1, as we reconstruct an Absorption Volume
		myMask = new DarkField3DTensorVolume(imgSizeX, imgSizeY, imgSizeZ, 1, getSpacing(), getOrigin());
		
		// Loop through all voxel elements
		for(int x = 0; x < imgSizeX; x++){
			for(int y = 0; y < imgSizeY; y++){
				for(int z = 0; z < imgSizeZ; z++){
					float val = reconAMP.getAtIndex(x, y, z, 0);
					// If value is in bounds set Mask to 1
					if(val > th_lower && val < th_higher){
						myMask.setAtIndex(x, y, z, 0, 1);
					} else{ // if value is out of bounds set Mask to 0
						myMask.setAtIndex(x, y, z, 0, 0);
					}
				} // END z
				} // END y
				} // END x
		
		return myMask;
	}

	
	/**
	 * Backprojects a sinogram into volume space
	 * A Parallel Beam is assumed. Also an easy trajectory without any deviations is assumed
	 * @param sino3D
	 * @return
	 */
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
