// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//
package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import weka.core.pmml.jaxbbindings.SetPredicate;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.data.numeric.Grid4D;
import edu.stanford.rsl.conrad.data.numeric.MultiChannelGrid3D;
import edu.stanford.rsl.conrad.numerics.SimpleVector;



public class DarkFieldGrid3DTensor extends Grid4D {
	
	// Defines dimension of volume box
	public int imgSizeX;
	public int imgSizeY;
	public int imgSizeZ;
	

	
	// CONSTRUCTOR for DarkFieldGrid3DTensor
	// width, height, depth gives the size of the volume
	// numChannels takes the number of vector elements to be saved in one voxel element
	
	public DarkFieldGrid3DTensor(int imgSizeX, int imgSizeY, int imgSizeZ, int numChannels){
		// Call constructor of MultiChannelGrid3D
		super(imgSizeX, imgSizeY, imgSizeZ, numChannels);
		this.imgSizeX=imgSizeX;
		this.imgSizeY=imgSizeY;
		this.imgSizeZ=imgSizeZ;
	
		
		
	}
	
	public void write3DTensorToImage(){
		
		
		
	}

			
	
	// Multiplies the whole grid with a given factor
	public void multiply(float factor){
		
		for(int x = 0; x <this.getSize()[0]; x++){
			for(int y = 0; y <this.getSize()[1]; y++){
				for(int z = 0; z <this.getSize()[0]; z++){
					multiplyAtIndexDarkfield(x,y,z,factor);
				} // End loop z
			} // End loop y
		} // End loop z
		
	}
	
	
	public double[] indexToPhysical(double x, double y, double z) {
		return new double[] { x * this.spacing[0] + this.origin[0],
				y * this.spacing[1] + this.origin[1],
				z * this.spacing[2] + this.origin[2]
		};
	}
	
	public int getNumberOfChannels(){
		// Return size of the last element of getSize()
		// This is the number of channels
		return this.getSize()[3];
	}
	
	
	// Multiplies the given grid tensor vector with a scalar value factor
	public void multiplyAtIndexDarkfield(int x, int y, int z, float factor){
		// Loop through all channels and multiply with factor
		for (int channel = 0; channel < this.getNumberOfChannels() ; channel++){
			// Multiply current value with factor
			float value = getAtIndex(x,y,z,channel);
			value = value*factor;
			setAtIndex(x,y,z,channel,value);
			
		}		
	}

	
	// Multiplies the entry at a given channel with a scalar value
	public void multiplyAtIndexDarkfield(int x, int y, int z, int channel, float val){
		setAtIndex(x, y,z, channel,getAtIndex(x, y, z, channel) *val);
	}
	
	
	
	// Store a channel value to given channel 
	public void setValueAtChannelN(int x, int y, int z, int channel, float value){
		setAtIndex(x,y,z,channel,value);
	}

	// Adds a given value to point at channel n 
	public void addAtIndexDarkfield(int x, int y, int z, int channel, float val){
		setAtIndex(x, y,z, channel,getAtIndex(x, y,z, channel) + val);
	}  

	
	// Store a whole tensor vector into grid point
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
	
	// Store a whole tensor vector into grid point
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
	
	
	public float[] getVectorAtIndex(int x, int y, int z) {
		
		float[] myVec = new float[getNumberOfChannels()];
		
		for (int channel = 0; channel < getNumberOfChannels(); channel++){
					myVec[channel] = getAtIndex(x, y, z, channel);
	}
		return myVec;
	}
	
	
	// Careful overrides ALL channels
	public void setAtIndex(int i, int j, int k, float val) {
		for(int channel = 0; channel < this.getNumberOfChannels(); channel++){
			setValueAtChannelN(i, j, k, channel, val);
		}
	}
	
	
	// Add a given tensor vector ontop of a grid point
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
	
	
	// Add a given tensor vector ontop of a grid point
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


	
	
	public void showComponents(){
		
		for(int channel = 0; channel < getNumberOfChannels(); channel++){
			String myTitle = "Volume at channel " + channel; 
			getSubGrid(channel).show(myTitle);
			
			
		}
	}
	
	
}



