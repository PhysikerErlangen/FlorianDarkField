// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//
package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.IJ;
import ij.ImagePlus;
import weka.core.pmml.jaxbbindings.SetPredicate;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.data.numeric.Grid4D;
import edu.stanford.rsl.conrad.data.numeric.MultiChannelGrid3D;
import edu.stanford.rsl.conrad.numerics.SimpleVector;


/**
 * 
 * 
 * @author Florian Schiffers
 *
 */
public class DarkFieldGrid3DTensor extends Grid4D {
	
	// Defines dimension of volume box
	public int imgSizeX;
	public int imgSizeY;
	public int imgSizeZ;
	
	
	/**
	 * @param imgSizeX
	 * @param imgSizeY
	 * @param imgSizeZ
	 * @param numChannels
	 */
	public DarkFieldGrid3DTensor(int imgSizeX, int imgSizeY, int imgSizeZ, int numChannels){
		// Call constructor of MultiChannelGrid3D
		super(imgSizeX, imgSizeY, imgSizeZ, numChannels);
		this.imgSizeX=imgSizeX;
		this.imgSizeY=imgSizeY;
		this.imgSizeZ=imgSizeZ;
	}
	
	
	/**
	 * write3DTensor writes the volume to a given file specified by filePath
	 * The name of the Image is not further defined, so default ""
	 * @param filePath - Path where Volume should be saved
	 */
	public void write3DTensorToImage(String filePath){
		write3DTensorToImage(filePath, "");
	}
	

	/**
	 * @param filePath - Path where Volume should be saved
	 * @param volumeName - Name of the image
	 */
	public void write3DTensorToImage(String filePath, String volumeName){
		
		ImagePlus myImage = edu.stanford.rsl.conrad.utils.ImageUtil.wrapGrid4D(this, volumeName);
		
		IJ.save(myImage,filePath);
		
	}


	/**
	 * Multiplies the whole grid with a given factor
	 * @param factor
	 */
	public void multiply(float factor){
		
		for(int x = 0; x <this.getSize()[0]; x++){
			for(int y = 0; y <this.getSize()[1]; y++){
				for(int z = 0; z <this.getSize()[0]; z++){
					multiplyAtIndexDarkfield(x,y,z,factor);
				} // End loop z
			} // End loop y
		} // End loop z
		
	}
	
	
	/**
	 * indexToPhysical calculates the world coordinate of a voxel element at (x,y,z)
	 * @param x
	 * @param y
	 * @param z
	 * @return Array of physical world coordinates of one voxel element
	 */
	public double[] indexToPhysical(double x, double y, double z) {
		return new double[] { x * this.spacing[0] + this.origin[0],
				y * this.spacing[1] + this.origin[1],
				z * this.spacing[2] + this.origin[2]
		};
	}
	
	/**
	 * @return Number of channels (e.g. scatterDirections)
	 */
	public int getNumberOfChannels(){
		// Return size of the last element of getSize()
		// This is the number of channels
		return this.getSize()[3];
	}
	
 
	/**
	 * Multiplies the complete given grid tensor vector with a scalar value factor
	 * @param x
	 * @param y
	 * @param z
	 * @param factor
	 */
	public void multiplyAtIndexDarkfield(int x, int y, int z, float factor){
		// Loop through all channels and multiply with factor
		for (int channel = 0; channel < this.getNumberOfChannels() ; channel++){
			// Multiply current value with factor
			float value = getAtIndex(x,y,z,channel);
			value = value*factor;
			setAtIndex(x,y,z,channel,value);
			
		}		
	}


	/**
	 * Multiplies the entry of a given channel with a scalar value
	 * @param x
	 * @param y
	 * @param z
	 * @param channel
	 * @param val
	 */
	public void multiplyAtIndexDarkfield(int x, int y, int z, int channel, float val){
		setAtIndex(x, y,z, channel,getAtIndex(x, y, z, channel) *val);
	}
	
 
	/**
	 * Stores a channel value to given channel
	 * @param x
	 * @param y
	 * @param z
	 * @param channel
	 * @param value
	 */
	public void setValueAtChannelN(int x, int y, int z, int channel, float value){
		setAtIndex(x,y,z,channel,value);
	}

	 
	/**
	 * Adds a given value to point at channel n
	 * @param x
	 * @param y
	 * @param z
	 * @param channel
	 * @param val
	 */
	public void addAtIndexDarkfield(int x, int y, int z, int channel, float val){
		setAtIndex(x, y,z, channel,getAtIndex(x, y,z, channel) + val);
	}  

	
	
	/**
	 * Stores a whole tensor vector into grid point
	 * @param x
	 * @param y
	 * @param z
	 * @param values
	 */
	public void setDarkFieldScatterTensor(int x, int y, int z, float[] values) {
		
		// Check if dimension of input and to be saved values are the same
		if( values.length != this.getNumberOfChannels() ){
		 throw new ArithmeticException("Dimension of input vector and vector to be set not equal.");
		}
		
		// Loop through all channels and set value
		for (int channel = 0; channel < values.length; channel++){
			setAtIndex(x,y,z,channel,values[channel]);
		}
	}
	

	/**
	 * stores a tensor vector into the given grid point at (x,y,z)
	 * @param x
	 * @param y
	 * @param z
	 * @param values
	 */
	public void setDarkFieldScatterTensor(int x, int y, int z, SimpleVector values) {
		
		// Check if dimension of input and to be saved values are the same
		if( values.getLen() != this.getNumberOfChannels() ){
		 throw new ArithmeticException("Dimension of input vector and vector to be set not equal.");
		}
		
		// Loop through all channels and set value
		for (int channel = 0; channel < values.getLen(); channel++){
			setAtIndex(x,y,z,channel,(float)values.getElement(channel));
		}
	}
	
	
	/**
	 * returns the tensor vector at a given voxel point
	 * @param x
	 * @param y
	 * @param z
	 * @return - vector
	 */
	public float[] getVectorAtIndex(int x, int y, int z) {
		
		float[] myVec = new float[getNumberOfChannels()];
		
		for (int channel = 0; channel < getNumberOfChannels(); channel++){
					myVec[channel] = getAtIndex(x, y, z, channel);
	}
		return myVec;
	}
	

	/**
	 * Set's the all values of the tensor vector to a specific value
	 * @param x
	 * @param y
	 * @param z
	 * @param val
	 */
	public void setAtIndex(int x, int y, int z, float val) {
		for(int channel = 0; channel < this.getNumberOfChannels(); channel++){
			setValueAtChannelN(x, y, z, channel, val);
		}
	}
	
	
	
	/**
	 * Add a given tensor vector ontop of a grid point
	 * @param x
	 * @param y
	 * @param z
	 * @param values
	 */
	public void subAtDarkFieldScatterTensor(int x, int y, int z, float[] values) {
		
		// Check if dimension of input and to be saved values are the same
		if( values.length != this.getNumberOfChannels() ){
		 throw new ArithmeticException("Dimension of input vector and vector to be set not equal.");
		}
			
		// Loop through all channels and add value
		for (int channel = 0; channel < values.length; channel++){
			setAtIndex(x,y,z,channel,getAtIndex(x, y,z, channel) - values[channel]);
		}
	}
	
	
	
	/**
	 * Adds a given tensor vector on top of a grid point
	 * @param x 
	 * @param y
	 * @param z
	 * @param values
	 */
	public void addAtDarkFieldScatterTensor(int x, int y, int z, float[] values) {
		
		// Check if dimension of input and to be saved values are the same
		if( values.length != this.getNumberOfChannels() ){
		 throw new ArithmeticException("Dimension of input vector and vector to be set not equal.");
		}
			
		// Loop through all channels and add value
		for (int channel = 0; channel < this.getNumberOfChannels(); channel++){
			float val = getAtIndex(x, y, z, channel) + values[channel];
			setAtIndex(x,y,z,channel,val);
		}
	}

	/**
	 * Displays each scatter direction component of the Grid4D as an own image
	 * All of them are Grid3D, which then can be displayed as a volume
	 * by the volume viewer.
	 * CAUTION: Auto contrast might be necessary
	 */
	public void showComponents(){
		
		for(int channel = 0; channel < getNumberOfChannels(); channel++){
			String myTitle = "Volume at channel " + channel; 
			getSubGrid(channel).show(myTitle);
		}
	}
	
/**
 * Displays a specific scatter direction component (this is a Grid3D)
 * @param channel
 */
public void showComponents(int channel){
			String myTitle = "Volume at channel " + channel; 
			getSubGrid(channel).show(myTitle);
	}
	
	
}



