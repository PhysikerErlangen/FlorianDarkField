//This code was developed in a collaboration with ECAP, Erlangen, Germany.
//This part of the code is not to be published under GPL before Oct 31st 2017.
//author@ Florian Schiffers July 1st, 2015
//
package edu.stanford.rsl.science.darkfield.FlorianDarkField;
import java.io.File;
import java.util.ArrayList;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.data.numeric.Grid4D;
import edu.stanford.rsl.conrad.geometry.General;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.ImageUtil;


public class DarkField3DTensorVolume extends DarkFieldGrid3DTensor{

	
	
	
	@SuppressWarnings("unused")
	private String title;
	
	
	/**
	 * @param imgSizeX - [px]
	 * @param imgSizeY - [px]
	 * @param imgSizeZ - [px]
	 * @param numChannels 
	 * @param spacing_world - [mm or other real world units]
	 * @param origin_world - [mm or other real world units]
	 */
	public DarkField3DTensorVolume(int imgSizeX,int imgSizeY, int imgSizeZ, int numChannels,double[] spacing_world, double[] origin_world){
		this(imgSizeX, imgSizeY, imgSizeZ, numChannels, spacing_world, origin_world, "");
	}
	
	/**
	 * @param imgSizeX - [px]
	 * @param imgSizeY - [px]
	 * @param imgSizeZ - [px]
	 * @param numChannels 
	 * @param spacing_world - [mm or other world units]
	 * @param origin_world - [mm or other world units]
	 * @param title
	 */
	public DarkField3DTensorVolume(int imgSizeX,int imgSizeY, int imgSizeZ, int numChannels,double[] spacing_world, double[] origin_world, String title) {
		
		// Call superconstructor of DarkFieldGrid3DTensor
		super(imgSizeX, imgSizeY, imgSizeZ,numChannels);
		
		// Set spacing of box
		setSpacing(spacing_world);
		// Set origin of 3D Image Box
		setOrigin(origin_world);
		// set Title
		setTitle(title);
	}

	
		
	public static void main(String[] args){
		
		String imagePath = "C:\\Users\\schiffers\\workspace\\MeasuredData\\DCI_volume.tif"; 
		 
		DarkField3DTensorVolume myVolume = readFromImagePlus(imagePath);
		new ImageJ();
		myVolume.show("Test Volume");
		myVolume.showComponent(0);
		
	}
	
	/**
	 * Reads the DarkField3DTensorVolume out of an ImagePlus
	 * @param imagePath
	 * @return
	 */
	public static DarkField3DTensorVolume readFromImagePlus(String imagePath){
		 ImagePlus imgVolume = IJ.openImage(imagePath);
		 return  readFromImagePlus(imgVolume);
	}
	
	
	/**
	 * @param imgVolume
	 * @return
	 */
	public static DarkField3DTensorVolume readFromImagePlus(ImagePlus imgVolume){
		
		// Open and read the image
		 
		// Make sure image is not null
		assert(imgVolume != null) : new Exception("Image could not be loaded.");
		
		int imgSizeX = imgVolume.getWidth();
		int imgSizeY = imgVolume.getHeight();
		int imgSizeZ = imgVolume.getNSlices();
		
		int numChannels = imgVolume.getNFrames();
		
		Calibration myCalib = imgVolume.getCalibration();
		
		double spacingX =  myCalib.pixelWidth;
		double spacingY =  myCalib.pixelHeight;
		double spacingZ =  myCalib.pixelDepth;
		
		double[] spacing = {spacingX,spacingY, spacingZ};
		
		double originX_pixel = myCalib.xOrigin;
		double originY_pixel = myCalib.yOrigin;
		double originZ_pixel = myCalib.zOrigin;
		
		double[] origin_world = {-originX_pixel*spacingX,-originY_pixel*spacingY,-originZ_pixel*spacingZ};
				
		DarkField3DTensorVolume myVolume = new DarkField3DTensorVolume(imgSizeX, imgSizeY, imgSizeZ, numChannels, spacing, origin_world);

		
		ImageStack myStack = imgVolume.getStack();
		
		
				for(int z = 0; z < imgSizeZ; z++){
					// Extracts the z Slice
					for(int channel = 0; channel < numChannels; channel++){
						// Calculates Stack Index
						int stackIndex = imgVolume.getStackIndex(0, z, channel);
						float[] sliceValues = ( float[]) myStack.getPixels(stackIndex);			
		
						for(int x = 0; x < imgSizeX; x++){
							for(int y = 0; y < imgSizeY; y++){
								int indexOnSlice = y*imgSizeX+x;
								float val = sliceValues[indexOnSlice];
								myVolume.setAtIndex(x, y, z, channel, val);
					} // END Y
				} // END X
			} // END CHANNEL
		} // END Z
		
		return myVolume;
	}
	
		
	
//	@Override
//	public Grid3D getChannel(int c){
//		Grid3D myGrid =  this.getSubGrid(c);
//		myGrid.setSpacing(getSpacing());
//		myGrid.setOrigin(getOrigin());
//		return myGrid;
//	}
//	
	
	public void setTitle(String t) {
		this.title = t;
	}
	
 
	/**
	 * masks a volume with a given mask. Careful, has to have same dimension
	 * @param mask
	 */
	public void maskWithVolume(DarkField3DTensorVolume mask){
		
		// If there's no mask just return and do nothing
		if(mask == null)
			return;
		
		for(int x = 0; x < imgSizeX; x++){
			for(int y = 0; y < imgSizeY; y++){
				for(int z = 0; z < imgSizeZ; z++){
					// The mask may only have ! one ! channel as it is as a mask
					if(mask.getAtIndex(x, y, z, 0) == 0f){
					this.setAtIndex(x, y, z, 0f);
					}
			}
		}
	}
	}		
	
	
	
	/**
	 * @param B
	 */
	public void sub( DarkField3DTensorVolume B){

		// Check for inconsistency (different dimensions)
		assert(getSize()[0]==B.getSize()[0]
				&&getSize()[1]==B.getSize()[1]
				&&getSize()[2]==B.getSize()[2])
				: new Exception("Dimension of data is wrong.");
		
		for(int x = 0; x <imgSizeX; x++){
			for(int y = 0; y <imgSizeY; y++){
				for(int z = 0; z <imgSizeZ; z++){
					
				this.subAtDarkFieldScatterTensor(x, y, z, B.getVectorAtIndex(x, y, z));
				
				} // End loop z
			} // End loop y
		} // End loop z
		
	}
	

	/**
	 * Return 
	 * @param A
	 * @param B
	 * @return result = A - B;
	 */
	public static DarkField3DTensorVolume sub(  DarkField3DTensorVolume A, DarkField3DTensorVolume B){

		// Check for inconsistency (different dimensions)
		assert(A.imgSizeX == B.imgSizeX
				&&A.imgSizeY == B.imgSizeY
				&&A.imgSizeZ == B.imgSizeZ
				): new Exception("Dimension of data is wrong.");
		
		DarkField3DTensorVolume sub = new DarkField3DTensorVolume(A.imgSizeX,A.imgSizeY, A.imgSizeX,
				A.getNumberOfChannels(), A.spacing, A.origin);
		
		for(int x = 0; x <A.imgSizeX; x++){
			for(int y = 0; y <A.imgSizeY; y++){
				for(int z = 0; z <A.imgSizeZ; z++){
					
					SimpleVector vecA = A.getSimpleVectorAtIndex(x, y, z);
					SimpleVector vecB = B.getSimpleVectorAtIndex(x, y, z);
					SimpleVector values = SimpleOperators.subtract(vecA, vecB);
					sub.setDarkFieldScatterTensor(x, y, z, values);
				} // End loop z
			} // End loop y
		} // End loop z
	
		return sub;
	}
	
	

	/**
	 * Class that calls the fiberDirectionWriterClass 
	 * @param pathFile
	 * @param fileName
	 */
	public void saveFiberOrientations(File pathFile, String fileName){
		DarkFieldReconPipeline.calculateFiberOrientations(this, DarkFieldScatterDirection.getScatterDirectionMatrix(this.getNumberOfChannels()), pathFile, fileName);
	}
	
	/**
	 * @author schiffers
	 * Different type of constraint published in paper of J. Vogel in
	 * Constrained X-Ray tensor tomography
	 */
	public static enum TensorConstraintType {
		/** HARD CONSTRAINT */
		HARD_CONSTRAINT,
		/** SOFT CONSTRAINT */ //TODO Implement
		SOFT_CONSTRAINT,
		/** NO CONSTRAINT AT ALL */ //TODO Implement
		NO_CONSTRAINT
	}
	
	
	public void enforceConstraint(TensorConstraintType type){
		if(type == TensorConstraintType.HARD_CONSTRAINT ){
			System.out.println("Enforce Hard Constraint on reconstructed scatter coefs: ");
			enforceHardConstraint();
		} else if (type == TensorConstraintType.SOFT_CONSTRAINT){
			// TODO NOT YET IMPLEMENTED
		} else{
			// TODO NOT YET IMPLEMENTED
		}
	} 
	
	private void enforceHardConstraint(){
		
		DarkFieldScatterExtractor myScatterExtractor = new DarkFieldScatterExtractor(this, DarkFieldScatterDirection.getScatterDirectionMatrix(this.getNumberOfChannels()));
	
		for(int x = 0; x <imgSizeX; x++){
			for(int y = 0; y <imgSizeY; y++){
			for(int z = 0; z <imgSizeZ; z++){
				SimpleVector constrainedCoef = myScatterExtractor.calcProjectedScatterCoefficients(x, y, z);
				this.setDarkFieldScatterTensor(x, y, z, constrainedCoef);
					
				}
			}
		}
		
				
	}

	
//	
//	@Override
//	public void show(){
//		show("");
//		
//	}
//	
//	@Override
//	public void show(String reconTitle){
//		
//		this.getMultichannelData().show();
//		
//		
//		for(int channel = 0; channel < this.getNumberOfChannels(); channel++){
//			
//			Grid3D curGrid = this.getChannel(channel);
//			
//			String title = reconTitle + "Channel Nr. "+channel;
//			ImagePlus myImg = ImageUtil.wrapGrid3D(curGrid, title);
//			
//			myImg.show();
//		
//			IJ.run(myImg, "Volume Viewer", "");
//			// curGrid.show();
//		}
//			
//			
//	}
	
}
