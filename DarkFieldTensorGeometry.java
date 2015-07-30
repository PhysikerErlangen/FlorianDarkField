// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;


import edu.stanford.rsl.conrad.geometry.Projection.CameraAxisDirection;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Point3D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.trajectories.Trajectory;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;



public class DarkFieldTensorGeometry {
	
	// maximum projection angle in index
	int maxTheta_index;   // [rad]
	// angle step size
	double deltaTheta; // [rad]
	
	// maximum Detector size
	int maxU_index;	  
	int maxV_index;	   
	
	// Step size projector
	double deltaU;	   // [mm]
	double deltaV;	   // [mm]
	//
	double maxU_world;
	double maxV_world;
	

	
	// Image Volume dimensions
	int imgSizeX;
	int imgSizeY;
	int imgSizeZ;
	
	// Image Volume dimensions
	double imgSizeX_world;
	double imgSizeY_world;
	double imgSizeZ_world;
	
	// Geometrical spacing of the volume [mm]
	double spacingX;
	double spacingY;
	double spacingZ;
	// 
	double originX;
	double originY;
	double originZ;

	// Magic offset TODO (I dont know what this is?!
	double offSetU_index;	   
	double offSetV_index;	   
	
	// DetectorOffset Given in Pixel
	double offSetU_world;	   
	double offSetV_world;
	
	// Number of scatter vectors
	int numScatterVectors;
		
	Configuration conf;
	
	Trajectory geo;
	
	// Permutation matrix, used to calculate ray direction for different trajectories
	String trajectoryFlag; // TODO initialize
	
	public DarkFieldTensorGeometry(Configuration conf, int numScatterVectors){
		
		this.conf = conf;
		
		this.geo = conf.getGeometry();
		
		maxU_index = geo.getDetectorWidth();
		maxV_index = geo.getDetectorHeight();
		
		deltaU = geo.getPixelDimensionX();
		deltaV = geo.getPixelDimensionY();
		
		maxU_world = maxU_index*deltaU;
		maxV_world = maxV_index*deltaV;
		
		imgSizeX = geo.getReconDimensionX();
		imgSizeY = geo.getReconDimensionY();
		imgSizeZ = geo.getReconDimensionZ();
		
		spacingX = geo.getVoxelSpacingX();
		spacingY = geo.getVoxelSpacingY();
		spacingZ = geo.getVoxelSpacingZ();
		
		
		imgSizeX_world = geo.getReconDimensionX()*spacingX;
		imgSizeY_world = geo.getReconDimensionY()*spacingY;
		imgSizeZ_world = geo.getReconDimensionZ()*spacingZ;
		
		
		// This is the actual origin of our Bounding Box
		originX = geo.getOriginX();
		originY = geo.getOriginY();
		originZ = geo.getOriginZ();
		
		// Calculate offset of detector in Pixel coordinates [px]
		offSetU_index = geo.getDetectorOffsetU(); //[px]
		offSetV_index = geo.getDetectorOffsetV(); //[px]
		
		// Calculate detector offset in world coordinates [mm]
		offSetU_world = offSetU_index*deltaU;	   
		offSetV_world = offSetV_index*deltaV;
		
		maxTheta_index = geo.getNumProjectionMatrices();
		
		double deltaThetaInDegree = geo.getAverageAngularIncrement();
		deltaTheta = Math.toRadians(deltaThetaInDegree);
			
		this.numScatterVectors = numScatterVectors;
	
		initRotMatrix();
		
	}

	// Checks if to SimpleVector are the same
	public boolean checkEquality(SimpleVector v1, SimpleVector v2){
		for (int k = 0; k < v1.getLen(); k++){
			if( ( v1.getElement(k) - v2.getElement(k) ) != 0){
				return false;
			}
		}
		return true;
	}
	
	
	public void initRotMatrix(){
		
		SimpleVector axis001 = new SimpleVector(0,0,1);
		SimpleVector axis010 = new SimpleVector(0,1,0);
		
		SimpleVector rotAxis = geo.getRotationAxis();
		
		if (checkEquality(axis001,rotAxis) ){
			trajectoryFlag = "001";
		}
		else if(checkEquality(axis010, rotAxis)){
			trajectoryFlag = "010";			
		} else{
			
		}
		
	}
	
	// Calculates the detector column in pixel coordinates  
	public double calcU_index(double uWorld){
		double curU_index = uWorld/deltaU - offSetU_index + maxU_index/2.0; 
		return curU_index;
	}
	
	public double calcV_index(double vWorld){
		double curV_index = vWorld/deltaV - offSetV_index + maxV_index/2.0; 
		return curV_index;
		
	}
	
	public double calculateDetectorCoordinate(int curU_index){
		// Calculate distance from camera center and include possible offset
		double s = deltaU * curU_index + this.offSetU_world - maxU_world/2.0;
		return s;
	}
	
	
	
	// Calculates the parallel projection of a "voxel" coordinate
	// onto the detector column in world coordinates
	public double calculateDetectorRow(int z_index){
		double curHeight = z_index*spacingZ + originZ;
		double curV = curHeight/deltaV - offSetV_index + maxV_index/2.0;
		return curV;
	}

	public double calculateHeight(double curV){
		// Calculate distance from camera center and include a possible offset
		// double curHeight = deltaV  * curV + offSetV_world - maxV/2.0;
		double curHeight = deltaV  * curV + offSetV_world - maxV_world/2.0;
		
		return curHeight;
		
	}
	
	public double[] getOrigin(){
		double[] origin = {geo.getOriginX(),geo.getOriginY(),geo.getOriginZ() }; 
		return origin;
		
	}
	
	public double[] getSpacing(){
		double[] origin = {geo.getVoxelSpacingX(),geo.getVoxelSpacingY(),geo.getVoxelSpacingZ()}; 
		return origin;
	}
	
	// checks if the point that should be interpolated lies inside of the bounding box.
	// if its to far from it, don't consider it and do not interpolate therefore!
	public boolean checkIfPointIsInBox(double x_ind, double y_ind, double z_ind){
		
		if (x_ind + 1 > imgSizeX || y_ind + 1 >= imgSizeY || z_ind + 1 >= imgSizeZ
				|| x_ind < 0 || y_ind < 0 || z_ind < 0){
		return false;
		
		}else{
			return true;
		}
		
	}
	
	
	public void calculateRayDirection(double curTheta){
		
		
		
	}
	
	
	
	// THIS is a really ugly method, but seems to be the fastest (but dirty) implementation
	// With this we implement two trajectories!
	public PointND calculateRotatedVector(double u_worldX, double u_worldY, double v_worldZ){

		if(trajectoryFlag.equals("001")){
		// Calculate Point when Rotation Axis is standard y - axis
		return new PointND(u_worldX,u_worldY,v_worldZ);
	}	else if(trajectoryFlag.equals("010")){
		return new PointND(-v_worldZ,u_worldY,u_worldX);
	}else{
		return null;
	}
}
	
	
	
	// Calculate orthogonal projection onto arbitrary plane by formula given by
	// https://de.wikipedia.org/wiki/Orthogonalprojektion
	// Returns: image coordinates in world coordinates (need to be transformed to image coordinates later)
	public static SimpleVector calcDetectorCoordinates(SimpleVector x, SimpleVector uVec, SimpleVector vVec){
		// under the assumption that u and v is orthogonal
		double u = SimpleOperators.multiplyInnerProd(x,uVec);
		double v = SimpleOperators.multiplyInnerProd(x,vVec);
				
		SimpleVector res = new SimpleVector(u,v);
		return res; 
		
	}

	
	
	// Calculate orthogonal projection onto arbitrary plane by formula given by
	// https://de.wikipedia.org/wiki/Orthogonalprojektion
	public static SimpleVector calcOrthogonalProjection(SimpleVector x, SimpleVector uVec, SimpleVector vVec){
		
		double inner1 = SimpleOperators.multiplyInnerProd(x,uVec);
		double inner2 = SimpleOperators.multiplyInnerProd(x,vVec);
		
		uVec.multiplyBy(inner1);
		vVec.multiplyBy(inner2);
		
		SimpleVector res = SimpleOperators.add(uVec,vVec);
		return res; 
		
	}
	
	
}
