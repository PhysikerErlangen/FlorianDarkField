// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;



import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;


public class ParallelDarkFieldProjector3DTensor extends  DarkFieldTensorGeometry {
	

	DarkFieldScatterWeightsCalculator scatterCoefficients;
	
	public ParallelDarkFieldProjector3DTensor(Configuration config, DarkFieldScatterWeightsCalculator scatterCoefficients){
	
		// Call super constructor
		super(config,scatterCoefficients.numScatterVectors);
		
		this.scatterCoefficients = scatterCoefficients;
		
	}
	

	

	// Project the Volume onto one projection
	public DarkField3DSinogram projectRayDriven(DarkField3DTensorVolume grid) {

		// Create sinogram to be reconstructed
		
		DarkField3DSinogram sino = new DarkField3DSinogram(maxU_index,maxV_index,maxTheta_index); //
		
		final double samplingRate = 3.d; // # of samples per pixel
				
		// Set spacing of sinogram

		sino.setSpacing(deltaU,deltaV, deltaTheta);

		
		// set up image bounding box in World Coordinates
		Translation trans = new Translation(
				// Translate the box by half the measuring volume
				// NumberPixelX * spacingX / 2
				originX,
				// NumberPixelY * spacingY / 2
				originY,
				// NumberPixelZ * spacingZ / 2
				originZ
			);
		
		// Inverse the transformation 
		Translation vectorToGridOrigin = trans.inverse();

		// One Box of the whole grid centered at 0
//		Box b = new Box((grid.getSize()[0] * grid.getSpacing()[0]), (grid.getSize()[1] * grid.getSpacing()[1]), (grid.getSize()[2] * grid.getSpacing()[2]));
		Box b = new Box(imgSizeX_world,imgSizeY_world,imgSizeZ_world);
	
		// Translate box to the given translation vector in coordinates
		b.applyTransform(trans);

		
		// Loop over all angles for the projection
		for(int curTheta=0; curTheta<maxTheta_index; curTheta++){
		
			// compute current angle theta [rad] and angular functions.
			double theta = deltaTheta * curTheta;
			// angular functions, precalculate for efficiency reasons
		
			double cosTheta = Math.cos(theta);
			double sinTheta = Math.sin(theta);
			
			
			// Console ouput: Uncomment if you debug in the DarkFieldProjector
			System.out.println("Cur Proj: " +curTheta +"/" +maxTheta_index 
					+ " (" +((10000*curTheta/maxTheta_index)/100.0) +"% done.)");
			

			// Go through detector line
			for (int curU = 0; curU < maxU_index; ++curU) {

				// Calculate distance from camera center and include possible offset
				double s = calcU_world(curU);
				
				
				// Go through all slices
				for( int curV = 0; curV < this.maxV_index; curV ++){
				
					double curHeight = curV_world(curV);
									
					// compute two points on the line through s and theta
				
				// One point is on the detector
				PointND p1 = calculateRotatedVector(s * cosTheta, s * sinTheta, curHeight);
				// Second is p1 plus a perpendicular vector
				PointND p2 = calculateRotatedVector(-sinTheta + (s * cosTheta),
						(s * sinTheta) + cosTheta, curHeight);

				// set up line equation between those 2 points
				StraightLine line = new StraightLine(p1, p2);
				// compute intersections between bounding box and intersection line.
				ArrayList<PointND> points = b.intersect(line);

				// only if we have intersections
				if (2 != points.size()){//??
					if(points.size() == 0) {
						line.getDirection().multiplyBy(-1.d);
						points = b.intersect(line);
					}
					if(points.size() == 0)
						continue;
				}

				PointND start = points.get(0); // [mm]
				PointND end = points.get(1);   // [mm]

				// get the normalized increment
				SimpleVector increment = new SimpleVector(
						end.getAbstractVector());
				increment.subtract(start.getAbstractVector());
				double distance = increment.normL2();
				// Check if both intersection are too close to each other (e.g. edge)
				double eps = 0.0001;
				if(distance < eps){
					continue;
				}
				
				increment.divideBy(distance * samplingRate);
				
				// Initialize the line integral
				float lineIntegral = 0.0f;
				
				// We need to shift into "Object" Coordinate system in order to 
				// calculate the indices used for interpolation later on
				start = vectorToGridOrigin.transform(start);

				// compute the integral along the line.
				for (double t = 0.0; t <=distance * samplingRate; ++t) {
					PointND current = new PointND(start);
					current.getAbstractVector().add(increment.multipliedBy(t));

					
					double x_ind = current.get(0) / spacingX;
					double y_ind = current.get(1) / spacingY;
					double z_ind = current.get(2) / spacingZ;
				
					if (checkIfPointIsInBox(x_ind,y_ind,z_ind) == false){
						continue;
					}

					for (int scatterIndex = 0; scatterIndex < this.numScatterVectors; scatterIndex++){
						// Get weight for current projection
						double scatterWeight = scatterCoefficients.getWeight(curTheta, scatterIndex);
						// Add to line integral
						
						
						lineIntegral += scatterWeight*InterpolationOperators.interpolateLinear(grid.getSubGrid(scatterIndex), x_ind, y_ind, z_ind);
					} // End scatter loop
				} // End ray tracing loop

				// normalize by the number of interpolation points
				lineIntegral /= samplingRate;
				
				
				
				
				
//				if(Float.isNaN(lineIntegral)){
//					System.out.println("We have an unwanted NaN here.");
//				}
//				
				
				// write integral value into the sinogram.
				sino.setAtIndex(curU, curV, curTheta, lineIntegral);
			
				} // End Height loop
			} // End Detector width loop
		} // End angle projection loop
		return sino;
}
	


//// SHOULD NOT WORK
// Project the Volume onto one projection
public DarkField3DSinogram projectPixelDriven(DarkField3DTensorVolume darkFieldVolume,DarkField3DTensorVolume ampMASK) {
	
	boolean debug = false;
	
	// Create sinogram to be reconstructed
	DarkField3DSinogram sino = new DarkField3DSinogram(this.maxU_index,this.maxV_index,this.maxTheta_index);
	
for( int curTheta = 0; curTheta < maxTheta_index; curTheta++){
	
	if(debug){
	// Debug Output: Uncomment if you need to debug the Backprojector
	System.out.println("Cur Proj: " +curTheta +"/" +maxTheta_index 
			+ " (" +((10000*curTheta/maxTheta_index)/100.0) +"% done.)");
	}
	
	// compute current angle theta [rad] and angular functions.
	double theta = deltaTheta * curTheta;
	// angular functions, precalculate for efficiency reasons

	double cosTheta = Math.cos(theta);
	double sinTheta = Math.sin(theta);
	
	SimpleVector dirU = calculateRotatedVector(cosTheta, sinTheta,0).getAbstractVector();
	SimpleVector dirV = calculateRotatedVector(0,0,1.0f).getAbstractVector();
	
	// Loop through all 3 dimensions
	for ( int x = 0; x < imgSizeX; x++){
		for(int y = 0; y < imgSizeY; y++){
			for(int z = 0; z < imgSizeZ; z++){
				
			if((ampMASK!=null)){
					if(ampMASK.getAtIndex(x, y, z, 0) == 0){
					continue;
					}
			}
				
			// compute world coordinate of current pixel
			 double[] w = darkFieldVolume.indexToPhysical(x, y,z);
				
			// Calculate current voxel element
			SimpleVector voxel = new SimpleVector(w[0],w[1],w[2]);				
			
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
					)
				continue; // Do nothing if projected point does not lie on detector
			
			float val = 0;
			
			for (int scatterChannel = 0; scatterChannel < this.numScatterVectors; scatterChannel++){
				// Get weight for current projection
				double scatterWeight = scatterCoefficients.getWeight(curTheta, scatterChannel);
				
				// getValue of current voxel element
				val += darkFieldVolume.getAtIndex(x,y,z,scatterChannel)*scatterWeight;
			}
			
			
			Grid2D projectedImage = sino.getSubGrid(curTheta);
			InterpolationOperators.addInterpolateLinear(projectedImage,curU_index,curV_index,val);
			
			} // END Z				
		}	// END Y
	} // END Z
} // END THETA
	return sino;
	}
}



//	
//	
//	// Calculates 
//	public int[] calcProjectionGeometry(int projIdx){
//		
//		// maxU
//		// maxV
//		
//		// maxTheta
//		
//		// num scatter Directions
//		
//		int[] test = {1,2};
//	
//		return test;
//		
//	}
//	
//	public double projectRayDriven(DarkField3DTensorVolume grid, int projIdx) {
//		
//		// Create a Grid2D with max detector size and maximum Theta angle for rotation
//		Grid2D detectorImage = new Grid2D(maxU_index,maxV_index);
//		detectorImage.setSpacing(deltaU, deltaV);
//		
//		
//		// set up image bounding box in World Coordinates		
//		Translation trans = new Translation(
//				// Translate the box by half the measuring volume
//				// NumberPixelX * spacingX / 2
//				-(grid.getSize()[0] * grid.getSpacing()[0])/2,
//				// NumberPixelY * spacingY / 2
//				-(grid.getSize()[1] * grid.getSpacing()[1])/2,
//				// NumberPixelZ * spacingZ / 2
//				-(grid.getSize()[2] * grid.getSpacing()[2])/2
//			);
//		
//		// Inverse the transformation 
//		Translation inverse = trans.inverse();
//
//		// One Box of the whole grid
//				Box boundingBox = new Box(imgSizeX*spacingX,imgSizeY*spacingY,imgSizeZ*spacingZ);
//					// Translate box to the origin
//				boundingBox.applyTransform(trans);
//		
//		
//				// TODO has to be implemented
//				int[] projGeo = calcProjectionGeometry(projIdx);
//				int thetaIdx = projGeo[0];
//				int coordU = projGeo[1];
//				int coordV = projGeo[2];
//		
//				
//
//				// compute two points on the line through s and theta
//				// We use PointND for Points in 3D space and SimpleVector for directions.
//				PointND p1 = new PointND(1, 0 , 0);
//				PointND p2 = new PointND(2, 0 , 0);
//				
//				// TODO Rotation of p1 and p2 with rotation matrices of projection
//				
//				// set up line equation
//				StraightLine line = new StraightLine(p1, p2);
//				// compute intersections between bounding box and intersection line.
//				ArrayList<PointND> points = boundingBox.intersect(line);
//
//				// only if 2 intersections have been found proceed
//				if (2 != points.size()){
//					
//					// TODO WHY should you do this?
//					if(points.size() == 0) {
//						// Change ray direction if there has been no intersection
//						line.getDirection().multiplyBy(-1.d);
//						points = boundingBox.intersect(line);
//					}
//					
//					// If there is no intersection between image box and ray
//					// Should not happen in general case, as mostly every ray is hitting the box
//					if(points.size() == 0)
//						return 0;
//					
//				}
//				
//				
//	
//				PointND start = points.get(0); // [mm]
//				PointND end = points.get(1);   // [mm]
//
//				
//				
//				// Get the increment of the ray normal vector for one step
//				final double samplingRate = 3.d; 
//				SimpleVector increment = new SimpleVector(
//						end.getAbstractVector());				
//				increment.subtract(start.getAbstractVector());
//				double distance = increment.normL2();
//				increment.divideBy(distance * samplingRate);
//				// End calculation of increment vector
//				
//				double lineIntegral = .0;
//				
//				// Transform start vector to start of bounding box 
//				start = inverse.transform(start);
//				
//				// compute the integral along the line.
//				
//				// Start with different scatterDireciton so scatter channel only needs to be load a few times
//				for (int scatterDir = 0; scatterDir < this.numScatterVectors; scatterDir ++){
//				
//					// Get scatterChannel of current signal
//					Grid3D scatterChannel = grid.getChannel(scatterDir);
//					
//					// Get scattering weight for current projection
//					double scatterWeight = scatterCoefficients.getWeight(projIdx, scatterDir);
//					
//					// Go through the ray with a well defined samplingRate and do bilinear Interpolation in every step
//					
//					for (double t = 0.0; t < distance * samplingRate; ++t) {
//					
//					PointND currentRayPoint = new PointND(start);
//					currentRayPoint.getAbstractVector().add(increment.multipliedBy(t));
//					
//					double x = currentRayPoint.get(0);
//					double y = currentRayPoint.get(1);
//					double z = currentRayPoint.get(2);
//
//					if (grid.getSize()[0] <= x + 1
//							|| grid.getSize()[1] <= y + 1
//							|| x < 0 || y < 0)
//						continue;
//					
//						// Add to line integral
//						lineIntegral += scatterWeight*InterpolationOperators.interpolateLinear(scatterChannel, x, y, z);
//					}
//				}
//
//				// normalize by the number of interpolation points
//				lineIntegral /= samplingRate;
//				// write integral value into the sinogram.
//		
//		return lineIntegral;
//	}
//	
//
//	// Project the Volume onto one projection
//	public float projectRayDriven(DarkField3DTensorVolume grid, int u, int v, int theta){
//		
//			final double samplingRate = 3.d; // # of samples per pixel
//					
//			// Set spacing of sinogram
//			
//			// set up image bounding box in World Coordinates
//			Translation trans = new Translation(
//					// Translate the box by half the measuring volume
//					// NumberPixelX * spacingX / 2
//					-(grid.getSize()[0] * grid.getSpacing()[0])/2,
//					// NumberPixelY * spacingY / 2
//					-(grid.getSize()[1] * grid.getSpacing()[1])/2,
//					// NumberPixelZ * spacingZ / 2
//					-(grid.getSize()[2] * grid.getSpacing()[2])/2
//				);
//			
//			// Inverse the transformation 
//			Translation inverse = trans.inverse();
//
//			// One Box of the whole grid centered at 0
//			Box b = new Box((grid.getSize()[0] * grid.getSpacing()[0]), (grid.getSize()[1] * grid.getSpacing()[1]), (grid.getSize()[2] * grid.getSpacing()[2]));
//		
//			// Translate box to the given translation vector in coordinates
//			b.applyTransform(trans);
//
//
//
//			
//			// Loop over all angles for the projection
//			for(int curTheta=0; curTheta<maxThetaIndex; ++curTheta){
//
//			
//				// compute current angle theta [rad] and angular functions.
//				double theta = deltaTheta * curTheta;
//				// angular functions, precalculate for efficiency reasons
//				
//				double cosTheta = Math.cos(theta);
//				double sinTheta = Math.sin(theta);
//
//				// Go through detector line
//				for (int curU = 0; curU < maxU; ++curU) {
//
//					// Calculate distance from camera center and include possible offset
//					double s = deltaU * curU - maxU / 2.0 + this.offSetU;
//					
//					
//					// Go through all slices
//					for( int curV = 0; curV < this.maxV; curV ++){
//					
//						double curHeight = deltaV * curV - maxV/2.0 + this.offSetV;
//						
//					// compute two points on the line through s and theta
//					
//					// One point is on the detector
//					PointND p1 = new PointND(s * cosTheta, s * sinTheta, curHeight);
//					// Second is p1 plus a perpendicular vector
//					PointND p2 = new PointND(-sinTheta + (s * cosTheta),
//							(s * sinTheta) + cosTheta, curHeight);
//
//					// set up line equation between those 2 points
//					StraightLine line = new StraightLine(p1, p2);
//					// compute intersections between bounding box and intersection line.
//					ArrayList<PointND> points = b.intersect(line);
//
//					// only if we have intersections
//					if (2 != points.size()){//??
//						if(points.size() == 0) {
//							line.getDirection().multiplyBy(-1.d);
//							points = b.intersect(line);
//						}
//						if(points.size() == 0)
//							continue;
//					}
//
//					PointND start = points.get(0); // [mm]
//					PointND end = points.get(1);   // [mm]
//
//					// get the normalized increment
//					SimpleVector increment = new SimpleVector(
//							end.getAbstractVector());
//					increment.subtract(start.getAbstractVector());
//					double distance = increment.normL2();
//					increment.divideBy(distance * samplingRate);
//
//					float lineIntegral = 0;
//			
//					start = inverse.transform(start);
//
//					// compute the integral along the line.
//					for (double t = 0.0; t <=distance * samplingRate; ++t) {
//						PointND current = new PointND(start);
//						current.getAbstractVector().add(increment.multipliedBy(t));
//
//						double x = current.get(0) / grid.getSpacing()[0];
//						double y = current.get(1) / grid.getSpacing()[1];
//						double z = curHeight;
//					
//
//						if (grid.getSize()[0] <=x + 1
//								|| grid.getSize()[1] <= y + 1
//								|| x < 0 || y < 0)
//							continue;
//						
//						for (int scatterIndex = 0; scatterIndex < this.numScatterVectors; scatterIndex++){
//							// Get weight for current projection
//							double scatterWeight = scatterCoefficients.getWeight(curU,curV,curTheta, scatterIndex);
//							// Add to line integral
//							lineIntegral += scatterWeight*InterpolationOperators.interpolateLinear(grid.getChannel(scatterIndex), x, y, z);
//						} // End scatter loop
//					} // End ray tracing loop
//
//					// normalize by the number of interpolation points
//					lineIntegral /= samplingRate;
//					// write integral value into the sinogram.
//					sino.setAtIndex(curU, curV, curTheta, lineIntegral);
//				
//					} // End Height loop
//				
//				} // End Detector width loop
//				
//			} // End angle projection loop
//			return sino;
//	}
//
//	
//	

