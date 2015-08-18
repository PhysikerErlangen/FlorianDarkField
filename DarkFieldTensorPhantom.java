package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.util.ArrayList;

import weka.core.pmml.jaxbbindings.GaussianDistribution;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import edu.stanford.rsl.apps.gui.opengl.PointCloudViewer;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.conrad.utils.FileUtil;
import edu.stanford.rsl.conrad.utils.RegKeys;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.DarkFieldPhantom;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.ImageToGrid3D;
import edu.stanford.rsl.science.darkfield.iterative.GradientSolver;
import edu.stanford.rsl.science.darkfield.reconstruction.WriteToVectorField;
import edu.stanford.rsl.tutorial.basics.MHDImageLoader;
import edu.stanford.rsl.tutorial.basics.PointCloudMaker;



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

	final static SimpleVector fiberDirX = new SimpleVector(1f,0f,0f);
	final static SimpleVector fiberDirY = new SimpleVector(0f,1f,0f);
	
//	public static void main (String [] args) throws Exception{
//	
//		String fileNameConfig1 = "E:\\fschiffers\\MeasuredData\\Phantom2\\PhantomHalfLarge_unsymetric.xml";
//		// Load configuration wooden case
//		Configuration config = Configuration.loadConfiguration(fileNameConfig1);
//		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig1);
//		// Reset rotation axis for Config2
//		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
//		Configuration2.getGeometry().setRotationAxis(rotationAxis2);
//		
//		
//		DarkFieldTensorPhantom myPhantom = new DarkFieldTensorPhantom(config);
//
//		ImagePlus myImage = DarkField3DTensorVolume.wrapDarkFieldGrid3DTensorToImagePlus(myPhantom.phantom, "test");
//		
//		// Load ImageJ
//		new ImageJ();
//		
//		myImage.show();
//		
//		
////		new ImageJ();
////		
////		myPhantom.phantom.show();
////		
////		myPhantom.calculateDarkFieldProjection(config,Configuration2);
////		
////		myPhantom.getDarkFieldSinogram(0).show();
////		myPhantom.getDarkFieldSinogram(1).show();
////		
////		myPhantom.getDarkFieldSinogram(0).showSinogram();
////		myPhantom.getDarkFieldSinogram(1).showSinogram();
////		
////		
////		System.out.println("Phaton created");
//	}



	
	/**
	 * @param config
	 */
	public DarkFieldTensorPhantom(Configuration config){
		this(config,7);
	}
	

	
	/**
	 * @param config
	 * @param numScatterVectors
	 */
	public DarkFieldTensorPhantom(Configuration config, int numScatterVectors){

		// Call super constructor of TensorGeometry
		super(config,numScatterVectors);
		
		// Init Phantom Volume
		phantom = new DarkField3DTensorVolume(this.imgSizeX,this.imgSizeY, this.imgSizeZ, numScatterVectors,this.getSpacing(),this.getOrigin());
		
		// Init Mask! Caution! numScatterVectors equals 1 as we have an "absorption" volume
		phantomMask = new DarkField3DTensorVolume(this.imgSizeX,this.imgSizeY, this.imgSizeZ,1,this.getSpacing(),this.getOrigin());
		
		sinogramList = new ArrayList<DarkField3DSinogram>();
		
		scatterDirections = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		calcPhantom();
	}
	
	public static SimpleVector calculateEllipsoidFromFiberOrientation(SimpleVector fiberDir, SimpleMatrix scatterDirections){
	
		SimpleVector eigenValues = new SimpleVector(7,7,1);
		SimpleMatrix eigenVectors = null;
		
				
		if(SimpleOperators.equalElementWise(fiberDir, fiberDirX,0.01)){
			// I define the columns of the SimpleMatrix to be the eigenVectors
			double[][] myEigenVectors ={ {0,0,1},{1,0,0},{0,1,0}};
			eigenVectors = new SimpleMatrix(myEigenVectors);
			
		} else if(SimpleOperators.equalElementWise(fiberDir, fiberDirY,0.01)){
			// I define the columns of the SimpleMatrix to be the eigenVectors
			double[][] myEigenVectors ={ {1,0,0},{0,0,1},{0,1,0}};
			eigenVectors = new SimpleMatrix(myEigenVectors);
		}
		
		SimpleVector scatterCoef = new SimpleVector(1f,1f,1f,1f,1f,1f,1f); 
		
		DarkFieldEllipsoid myEllipsoid = new DarkFieldEllipsoid(scatterDirections, scatterCoef, eigenValues, eigenVectors);
		
		SimpleVector projScatterCoef = myEllipsoid.calculateSquaredProjectedCoefficients();
		
		return projScatterCoef;
	}
	
	
	/**
	 * 
	 */
	public void calcPhantom(){
		
		SimpleVector fiberDir1 = new SimpleVector(1f,0f,0f);
		SimpleVector fiberDir2 = new SimpleVector(0f,1f,0f);
		
		SimpleVector scatterCoefDir1 = calculateEllipsoidFromFiberOrientation(fiberDir1, scatterDirections);
		SimpleVector scatterCoefDir2 = calculateEllipsoidFromFiberOrientation(fiberDir2, scatterDirections);
		
		int aX = (int)( 0.3*imgSizeX);
		int bX = (int) (0.7*imgSizeX);
		
		int aY = (int)( 0.4*imgSizeY);
		int bY = (int) (0.6*imgSizeY);
		
		int l1 = (int)( 0.1*imgSizeZ);
		int l2 = (int) (0.3*imgSizeZ);
		int l3 = (int) (0.5*imgSizeZ);
		int l4 = (int) (0.7*imgSizeZ);
		int l5 = (int) (0.9*imgSizeZ);
		
		
		
		// System.out.println("a:" + aX + "b: " + bX);
		
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
		
		DarkField3DSinogram sino1 = projector1.projectPixelDriven(phantom);
		DarkField3DSinogram sino2 = projector2.projectPixelDriven(phantom);
		
		
		
		sinogramList.add(sino1);
		sinogramList.add(sino2);
		
	}
	
}