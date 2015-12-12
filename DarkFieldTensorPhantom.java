package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.util.ArrayList;

import ij.ImagePlus;
import edu.stanford.rsl.conrad.geometry.General;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;



public class DarkFieldTensorPhantom  extends  DarkFieldTensorGeometry  {

	
	private ArrayList<DarkField3DSinogram> sinogramList;
	private ArrayList<ImagePlus> projectionList;
	
	private DarkField3DTensorVolume phantomMask;
	
	/**
	 * @return the phantomMask
	 */
	public DarkField3DTensorVolume getPhantomMask() {
		return phantomMask;
	}



	/**
	 * @param phantomMask the phantomMask to set
	 */
	public void setPhantomMask(DarkField3DTensorVolume phantomMask) {
		this.phantomMask = phantomMask;
	}

	private DarkField3DTensorVolume phantom;
	
	/**
	 * @return the phantom
	 */
	public DarkField3DTensorVolume getPhantom() {
		return phantom;
	}

	/**
	 *  SCATTER DIRECTION MATRIX
	 */
	private SimpleMatrix scatterDirections;
	
	/**
	 * @return the scatterDirections
	 */
	public SimpleMatrix getScatterDirections() {
		return scatterDirections;
	}

	final static SimpleVector fiberDir1 = new SimpleVector(1f,0f,0f).normalizedL2();
	final static SimpleVector fiberDir2 = new SimpleVector(0f,1f,0f).normalizedL2();
	final static SimpleVector fiberDir3 = new SimpleVector(1f,1f,0f).normalizedL2();
	final static SimpleVector fiberDir4 = new SimpleVector(-1,1f,0f).normalizedL2();
	
	
	public static enum PhantomType {
		/**
		 * TODO DESCRIPTION
		 */
		WOODEN_BLOCK_PHANTOM,
		/**
		 * TODO DESCRIPTION 
		 */
		CURL_VECTOR_FIELD_PHANTOM,
		/**
		 * TODO DESCRIPTION 
		 */
		PATCHY_PHANTOM,
		/**
		 * TODO DESCRIPTION 
		 */
		WHATEVER_PHANTOM,
		/**
		 * TODO DESCRIPTION
		 */
		MY_CRAZY_PHANTOM2
	}
	
	
	
	/**
	 * CONSTRUCTORS
	 */
	
	
	/**
	 * @param config
	 * @param numScatterVectors
	 */
	public DarkFieldTensorPhantom(Configuration config, int numScatterVectors,PhantomType myPhantomType){

		// Call super constructor of TensorGeometry
		super(config,numScatterVectors);
		
		// Init Phantom Volume
		phantom = new DarkField3DTensorVolume(this.imgSizeX,this.imgSizeY, this.imgSizeZ, numScatterVectors,this.getSpacing(),this.getOrigin());
		
		// Init Mask! Caution! numScatterVectors equals 1 as we have an "absorption" volume
		phantomMask = new DarkField3DTensorVolume(this.imgSizeX,this.imgSizeY, this.imgSizeZ,1,this.getSpacing(),this.getOrigin());
		
		sinogramList = new ArrayList<DarkField3DSinogram>();
		
		scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		// calcWoodenBlockPhantom();
		calcPhantom(myPhantomType);
	}
	

	/**
	 * @param myPhantomType
	 */
	private void calcPhantom(PhantomType myPhantomType){
		
		switch(myPhantomType){
		case CURL_VECTOR_FIELD_PHANTOM:
			calcCurlyPhantom();
			break;
		case PATCHY_PHANTOM:
			calcPatchyPhantom();
			break;
		case WOODEN_BLOCK_PHANTOM:
			calcWoodenBlockPhantom();
			break;
		default:
		calcWoodenBlockPhantom();
		break;
		}
		
	}
	
	/**
	 * @param fiberDir
	 * @param scatterDirections
	 * @return
	 */
	public static SimpleVector calculateEllipsoidFromFiberOrientation(SimpleVector fiberDir, SimpleMatrix scatterDirections){

		// First normalize the vector in case this is not done yet
		fiberDir = fiberDir.normalizedL2();
		
		if(Double.isNaN(fiberDir.getElement(0))){
			fiberDir = new SimpleVector(0,0,0);
		}
		
		// helper is vector linear independent to fiberDir
		SimpleVector helper;
		// to ensure linear independence check first
		if(SimpleOperators.equalElementWise(fiberDir.normalizedL2().absoluted(),fiberDir1,0.0001)){
			 helper = fiberDir2;
		} else{
			helper = fiberDir1;
		}
		SimpleVector eigenVec1 = General.crossProduct(helper, fiberDir).normalizedL2();
		SimpleVector eigenVec2 = General.crossProduct(eigenVec1, fiberDir).normalizedL2(); 
		
		SimpleMatrix eigenVectors = new SimpleMatrix(3,3);
		
		eigenVectors.setColValue(0, eigenVec1);
		eigenVectors.setColValue(1, eigenVec2);
		eigenVectors.setColValue(2, fiberDir);
		
		int numScatterVectors = scatterDirections.getCols();
		
		SimpleVector scatterCoef = new SimpleVector(numScatterVectors);
		for(int channel = 0; channel < numScatterVectors; channel++){
			scatterCoef.setElementValue(channel, 1);
		}
		double ratio = 7;
		double lambda3 = 1;
		double lambda2 = ratio*lambda3;
		double lambda1 = ratio*lambda3;
		
		SimpleVector eigenValues = new SimpleVector(lambda1,lambda2,lambda3);
		
		DarkFieldEllipsoid myEllipsoid = new DarkFieldEllipsoid(scatterDirections, scatterCoef, eigenValues, eigenVectors);
		
		SimpleVector projScatterCoef = myEllipsoid.calculateSquaredProjectedCoefficients();
		
		return projScatterCoef;
	}
	
	
	/**
	 * @param imgSizeX
	 * @param imgSizeY
	 * @return
	 */
	private SimpleVector[][] calcCurlLayer(int imgSizeX, int imgSizeY){
		
		SimpleVector eX = new SimpleVector(1,0,0);
		SimpleVector eY = new SimpleVector(0,1,0);
		
		SimpleVector[][] myVectorField = new SimpleVector[imgSizeX][imgSizeY];
		
		for(int x = 0; x < imgSizeX; x++){
			for(int y = 0; y < imgSizeY; y++){
				/**
				 * Calculate vector field
				 */
				SimpleVector vec = SimpleOperators.subtract(eX.multipliedBy(y-imgSizeY/2.0),eY.multipliedBy(x-imgSizeX/2.0));
				
				vec.normalizeL2();
				// Vec is 0-Vector normalization will result in NaN, so we have to check for that
				if(Double.isNaN(vec.getElement(0))){
					vec = new SimpleVector(1,0,0);
				}
				
				/*
				 * Calculate scatter coefficients
				 */
				SimpleVector scatterCoef = calculateEllipsoidFromFiberOrientation(vec, scatterDirections);
				
				if(Double.isNaN(scatterCoef.getElement(0))){
					// TEst
					System.out.println("NaN - Entry. Something went wrong! Better check this!");
				}
				myVectorField[x][y] = scatterCoef;
				
			}
		}
		
		return myVectorField;
	}
	
	
	/**
	 * Calculates a Dark Field Phantom containg a curly layer
	 */
	public void calcCurlyPhantom(){
	
		SimpleVector[][] myCurlVector = calcCurlLayer(imgSizeX, imgSizeY);
		
		int aX = (int)( 0.3*imgSizeX);
		int bX = (int) (0.7*imgSizeX);
		
		int aY = (int)( 0.4*imgSizeY);
		int bY = (int) (0.6*imgSizeY);
		
		int l1 = (int)( 0.1*imgSizeZ);
		int l2 = (int) (0.3*imgSizeZ);
		int l3 = (int) (0.5*imgSizeZ);
		int l4 = (int) (0.7*imgSizeZ);
		int l5 = (int) (0.9*imgSizeZ);
		
		for(int x = aX; x < bX; x ++){
			for(int y = aY; y < bY; y ++){
				for(int z = l1; z <= l5; z ++){
					// Creates a phantom of 3 Layers, each Layer contains voxel with a constant scatter direction!
					if(z<l2){
						phantom.setDarkFieldScatterTensor(x, y, z, myCurlVector[x][y]);	
					} else if(z<l3){
						phantom.setDarkFieldScatterTensor(x, y, z, myCurlVector[x][y]);
					} else if(z<l4){
						phantom.setDarkFieldScatterTensor(x, y, z, myCurlVector[x][y]);
					} else if(z<l5){
						phantom.setDarkFieldScatterTensor(x, y, z, myCurlVector[x][y]);
					}  
					phantomMask.setDarkFieldScatterTensor(x, y, z, new SimpleVector(1f));
				} // END X 
		} // END Y
	} // END Z
	
		
	}
	
	/**
	 * 
	 */
	public void calcWoodenBlockPhantom(){
		
		SimpleVector scatterCoefDir1 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,0,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir2 = calculateEllipsoidFromFiberOrientation(new SimpleVector(0,1,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir3 = calculateEllipsoidFromFiberOrientation(new SimpleVector(0,0,1).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir4 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,1,1).normalizedL2(), scatterDirections);
		
		int aX = (int)( 0.1*imgSizeX);
		int bX = (int) (0.9*imgSizeX);
		
		int aY = (int)( 0.1*imgSizeY);
		int bY = (int) (0.9*imgSizeY);
		
		int l1 = (int)( 0.1*imgSizeZ);
		int l2 = (int) (0.3*imgSizeZ);
		int l3 = (int) (0.5*imgSizeZ);
		int l4 = (int) (0.7*imgSizeZ);
		int l5 = (int) (0.9*imgSizeZ);
		
		for(int x = aX; x < bX; x ++){
			for(int y = aY; y < bY; y ++){
				for(int z = l1; z <= l5; z ++){
					// Creates a phantom of 3 Layers, each Layer contains voxel with a constant scatter direction!
					if(z<l2){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir1);	
					} else if(z<l3){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir2);
					} else if(z<l4){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir3);
					} else if(z<l5){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir4);
					}  
					// Set Phantom Mask
					phantomMask.setDarkFieldScatterTensor(x, y, z, new SimpleVector(1f));
				} // END X 
		} // END Y
	} // END Z
	

} // End calcPhantom()
	
	
	private SimpleVector[][][] calcPathyDirections(){
		
		SimpleVector[][][] myDirMatrix = new SimpleVector[imgSizeX][imgSizeY][imgSizeZ];
		
//		for(int x = 0; x < 4; x++){
//			for(int y = 0; y < imgSizeY; y++){
//				for(int z = 0; z < imgSizeZ; z++){
//					double a =  
//					double b = 
//					double c = 
//					SimpleVector dir = new SimpleVector(a,b,c);
//						
//					}
//				}
//			}
		
		return myDirMatrix;
		
	}
	
	public void calcPatchyPhantom(){
		
		SimpleVector scatterCoefDir1 = calculateEllipsoidFromFiberOrientation(new SimpleVector(18,0,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir2 = calculateEllipsoidFromFiberOrientation(new SimpleVector(0,1,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir3 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,1,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir4 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,-1,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir5 = calculateEllipsoidFromFiberOrientation(new SimpleVector(2,0,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir6 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,0,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir7 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,0,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir8 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,0,0).normalizedL2(), scatterDirections);
		SimpleVector scatterCoefDir9 = calculateEllipsoidFromFiberOrientation(new SimpleVector(1,0,0).normalizedL2(), scatterDirections);
		
		int aX = (int)( 0.3*imgSizeX);
		int bX = (int) (0.7*imgSizeX);
		
		int aY = (int)( 0.4*imgSizeY);
		int bY = (int) (0.6*imgSizeY);
		
		int l1 = (int)( 0.1*imgSizeZ);
		int l2 = (int) (0.3*imgSizeZ);
		int l3 = (int) (0.5*imgSizeZ);
		int l4 = (int) (0.7*imgSizeZ);
		int l5 = (int) (0.9*imgSizeZ);
		
		for(int x = aX; x < bX; x ++){
			for(int y = aY; y < bY; y ++){
				for(int z = l1; z <= l5; z ++){
					// Creates a phantom of 3 Layers, each Layer contains voxel with a constant scatter direction!
					if(z<l2){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir1);	
					} else if(z<l3){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir2);
					} else if(z<l4){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir1);
					} else if(z<l5){
						phantom.setDarkFieldScatterTensor(x, y, z, scatterCoefDir2);
					}  
					// Set Phantom Mask
					phantomMask.setDarkFieldScatterTensor(x, y, z, new SimpleVector(1f));
				} // END X 
		} // END Y
	} // END Z
	

} // End calcPhantom()
	
	
	
	/**
	 * @param index
	 * @return
	 */
	public DarkField3DSinogram getDarkFieldSinogram(TrajectoryType type){
		if(type == TrajectoryType.HORIZONTAL){
			return sinogramList.get(0);
		} else if(type == TrajectoryType.VERTICAL){
			return sinogramList.get(1);
		} else {
			return null;
		}

	}
	
	/**
	 * @param config1
	 * @param config2
	 */
	public void calculateDarkFieldProjection(Configuration config1, Configuration config2){
		
		DarkFieldScatterWeightsCalculator scatterCoef1 = new DarkFieldScatterWeightsCalculator(config1,phantom.getNumberOfChannels());
		DarkFieldScatterWeightsCalculator scatterCoef2 = new DarkFieldScatterWeightsCalculator(config2,phantom.getNumberOfChannels());
		
		ParallelDarkFieldProjector3DTensor projector1 = new
				ParallelDarkFieldProjector3DTensor(config1, scatterCoef1);
		
		ParallelDarkFieldProjector3DTensor projector2 = new
				ParallelDarkFieldProjector3DTensor(config2, scatterCoef2);
		
		//DarkField3DSinogram sino1 = projector1.projectRayDriven(phantom);
		//DarkField3DSinogram sino2 = projector2.projectRayDriven(phantom);
		
		DarkField3DSinogram sino1 = projector1.projectPixelDriven(phantom,phantomMask);
		DarkField3DSinogram sino2 = projector2.projectPixelDriven(phantom,phantomMask);
		
		
		
		sinogramList.add(sino1);
		sinogramList.add(sino2);
		
	}
	
}