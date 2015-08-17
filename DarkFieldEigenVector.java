package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.data.numeric.Grid4D;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

/**
 * Defines three dimensional vector field
 *
 */
public class DarkFieldEigenVector extends Grid4D {


	// Defines dimension of volume box
	protected int imgSizeX;
	protected int imgSizeY;
	protected int imgSizeZ;
	
	/**
	 * @param imgSizeX
	 * @param imgSizeY
	 * @param imgSizeZ
	 * @param spacing_world
	 * @param origin_world
	 */
	public DarkFieldEigenVector(int imgSizeX, int imgSizeY,int imgSizeZ, double[] spacing_world, double[] origin_world){
		// Dimension of 4!
		// 
		super(imgSizeX, imgSizeY, imgSizeZ, 3);
		
		this.imgSizeX = imgSizeX;
		this.imgSizeY = imgSizeY;
		this.imgSizeZ = imgSizeZ;
		
		// Set spacing of box
		setSpacing(spacing_world);
		// Set origin of 3D Image Box
		setOrigin(origin_world);
		
	}
	
	/**
	 * @param x
	 * @param y
	 * @param z
	 * @param vec - EigenValue is encoded in length of vector
	 */
	public void setVector(int x, int y, int z, SimpleVector vec){
	
		// Check for inconsistency (different dimensions)
		assert(vec.getLen()==3): new Exception("Dimension of data vector has to be 3.");
		// Loop through every coordinate
		for(int i = 0; i < 3; i++){
			super.setAtIndex(x, y, z, i, (float) vec.getElement(i));
		}
	}
	
	/**
	 * returns the world vector at a given voxel point
	 * @param x
	 * @param y
	 * @param z
	 * @return - vector
	 */
	public SimpleVector getSimpleVectorAtIndex(int x, int y, int z) {

		SimpleVector myVec = new SimpleVector(3);
		
		for (int coord = 0; coord < 3; coord++){ // 3 world coordinates x,y and z
					myVec.setElementValue(coord, getAtIndex(x, y, z, coord));
	}
		return myVec;
	}
	
	
	
	
	
}
