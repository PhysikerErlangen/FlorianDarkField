// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//


package edu.stanford.rsl.science.darkfield.FlorianDarkField;


import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;









import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Transform;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;


public class ParallelDarkFieldBackprojector3DTensor extends  DarkFieldTensorGeometry{
	
	// Sampling rate of line integral
	final double samplingRate = 3.d;
	
	DarkFieldScatterCoef scatterCoefficients;
	
	
	
	/**
	 * Sampling of projections is defined in the constructor.
	 * 
	 * @param imageSizeU
	 * @param imageSizeV
	 * @param pxSzUMM
	 * @param pxSzVMM
	 */
	public ParallelDarkFieldBackprojector3DTensor(Configuration config, DarkFieldScatterCoef scatterCoefficients) {
	
		// Open super operator of geometry class
		super(config,scatterCoefficients.numScatterVectors);
		
		this.scatterCoefficients = scatterCoefficients;
		
	}


	
	
	
	/**
	 * The ray driven solution.
	 * 
	 * @param sino
	 *            the sinogram
	 * @return the image
	 */
	
	public DarkField3DTensorVolume backprojectRayDriven(DarkField3DSinogram sino3D) {
				
		// Create a clear 3D Tensor volume image
		DarkField3DTensorVolume darkFieldVolume = new DarkField3DTensorVolume(this.imgSizeX,
				this.imgSizeY, this.imgSizeZ,numScatterVectors,getSpacing(),getOrigin());

		// Vector operation that translates by origin of the volume to be reconstructed
		Translation trans = new Translation(getOrigin());
		
		// calculate transInverse translation vector
		Transform transInverse = trans.inverse();

		// Create a bounding box to calculate the crossing of the ray tracing lines later on 
		Box boundingBox = new Box(imgSizeX_world,imgSizeY_world,imgSizeZ_world);
		
		// Shift the bounding box to origin of the volume to be reconstructed
		boundingBox.applyTransform(trans);

		// Iterate over all projection angles
		for (int curTheta = 0; curTheta < maxTheta_index; ++curTheta) {
		
			// precompute theta [rad] and angular functions.
			double theta = deltaTheta * curTheta;
			double cosTheta = Math.cos(theta);
			double sinTheta = Math.sin(theta);

			// Debug Output: Uncomment if you need to debug the Backprojector
			System.out.println("Cur BackProj: " +curTheta +"/" +maxTheta_index 
					+ " (" +((10000*curTheta/maxTheta_index)/100.0) +"% done.)");
			
			// Loop through all elements of U-Axis of detector
			for (int curU = 0; curU < this.maxU_index; curU++) {
			
				// Calculate distance from camera center and include possible offset
				double s = calculateDetectorCoordinate(curU);
				
				// Loop through all slices of detector
				for(int curV = 0; curV < this.maxV_index; curV++){
				
				// Calculate vertical distance from camera center and include a possible offset
				double curHeight = calculateHeight(curV);
					
				// compute two points on the line through s and theta
				// We use PointND for Points in 3D space and SimpleVector for
				// directions.
				// One point lies on the detector
				PointND p1 = calculateRotatedVector(s * cosTheta, s * sinTheta, curHeight);
				// Second point is just added to p1 but slight perpendicular to it
				PointND p2 = calculateRotatedVector(-sinTheta + (s * cosTheta),
						(s * sinTheta) + cosTheta, curHeight);
				
				// set up line equation
				StraightLine line = new StraightLine(p1, p2);
				
				// compute intersections between bounding box and intersection line.
				ArrayList<PointND> points = boundingBox.intersect(line);

				// only if we have intersections
				if (2 != points.size()) {
					if (points.size() == 0) {
						line.getDirection().multiplyBy(-1.d);
						points = boundingBox.intersect(line);
					}
					if (points.size() == 0)
						continue;
				}
				
				PointND start = points.get(0);  //[mm]
				PointND end = points.get(1);    //[mm]

				// calculate the normalized increment
				SimpleVector increment = new SimpleVector(end.getAbstractVector());
				// Substract start from end to get a "direction" vector
				increment.subtract(start.getAbstractVector());
				// calculate norm of direction vector and normalize with respect to sampling rate
				double distance = increment.normL2();

				// Check if both intersection are too close to each other (e.g. edge)
				if(distance < 0.001){
					continue;
				}

				increment.divideBy(distance * samplingRate);

				// Get darkfield Value out of sinogram
				float darkFieldValue = sino3D.getAtIndex(curU, curV, curTheta);
								
				//TODO WHY TO DO THE INVERSE TRANSFORMATION ?!
				start = transInverse.transform(start);

				// compute the integral along the line.
				
				for (double t = 0.0; t < distance * samplingRate; ++t) {
					// Calculate current point
					PointND current = new PointND(start);
					current.getAbstractVector().add(increment.multipliedBy(t));

					// Calculate current XYZ - Point in World Coordinates
					
					/* Increment vector is multiplied with the samplingRate
					 * To calculate real world value we have to divide with spacing
					 */
					
					double x_ind = current.get(0) / spacingX;
					double y_ind = current.get(1) / spacingY;
					double z_ind = current.get(2) / spacingZ;

					if (checkIfPointIsInBox(x_ind,y_ind,z_ind) == false){
						continue;
					}

					
					// Loop through all scatter vectors					
					for ( int scatterChannel = 0; scatterChannel < this.numScatterVectors; scatterChannel++){
					
						// Get weight for current projection
						double scatterWeight = scatterCoefficients.getWeight(curTheta, scatterChannel);
										
					
					// Interpolate point
					InterpolationOperators.addInterpolateLinear(darkFieldVolume.getSubGrid(scatterChannel),x_ind, y_ind, z_ind, (float)(darkFieldValue*scatterWeight));
					
					} // END SCATTER LOOP
				} // END RAY TRACING LOOP
				}  // END HEIGT LOOP
			} // END LINE DETECTOR LOOP 
		} // END ANGLE LOOP
		

		
		// TODO Normalization FACTOR (WHAT IS THIS USED FOR )
		
		float normalizationFactor = (float) (samplingRate * maxTheta_index / deltaU / Math.PI);
		
		NumericPointwiseOperators.divideBy(darkFieldVolume, normalizationFactor);
		
		return darkFieldVolume;
	}

	
	
	public	DarkField3DTensorVolume backprojectPixelDriven(DarkField3DSinogram sino3D) {
		
		boolean debug = true;
		
		// Create empty 3DDarkField Volume
		DarkField3DTensorVolume grid = new DarkField3DTensorVolume(imgSizeX,imgSizeY,imgSizeZ,numScatterVectors,getSpacing(),getOrigin());
		
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
					float darkFieldValue = InterpolationOperators.interpolateLinear(detectorImageAtTheta,curU_index,curV_index);
					
					float[] values = new float[numScatterVectors];
					
					for (int scatterChannel = 0; scatterChannel < numScatterVectors; scatterChannel++){
						
						// Works only for parallel beam, as weight is the same for every parallel ray!
						double scatterWeight = scatterCoefficients.getWeight(curTheta, scatterChannel);
						float val = (float)(scatterWeight*darkFieldValue);
						values[scatterChannel] = val;
						
					} 
					
					grid.addAtDarkFieldScatterTensor(x, y, z, values);
					
					
					
					// END LOOP DARK FIELD SCATTERER
					} // END LOOP Z
						
						
						
				} // END LOOP Y
					
			} // END LOOP X
		} // END LOOP ANGLES
		
		
		
		// TODO Normalization factors
		NumericPointwiseOperators.divideBy(grid, (float) (maxTheta_index / Math.PI));
		return grid;
	}
	
	
	

}