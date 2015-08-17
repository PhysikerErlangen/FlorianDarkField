package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldScatterExtractor {

	
	private DarkField3DTensorVolume darkFieldVolume;
	private SimpleMatrix scatterMatrix;
	
	private DarkFieldTensorClass tensorClass;
	
	private int imgSizeX;
	private int imgSizeY;
	private int imgSizeZ;
	
	

	/**
	 * @param darkFieldVolume - reference to the reconstructed dark field volume
	 * @param scatterMatrix - Contains the scatter directions in a 3 X NumScattervectors Matrix
	 */
	public DarkFieldScatterExtractor(DarkField3DTensorVolume darkFieldVolume, SimpleMatrix scatterMatrix){
		
		this.darkFieldVolume = darkFieldVolume;
		this.scatterMatrix = scatterMatrix;
		
		this.imgSizeX = darkFieldVolume.imgSizeX;
		this.imgSizeY = darkFieldVolume.imgSizeY;
		this.imgSizeZ = darkFieldVolume.imgSizeZ;
		
	}
	
	

	
	
	/**
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	public SimpleVector calcProjectedScatterCoefficients(int x, int y, int z){
		

		// Get scatter coefficient of reconstruction
		SimpleVector scatterCoef = darkFieldVolume.getSimpleVectorAtIndex(x, y, z);

		// Initializes the PCA object
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterMatrix,scatterCoef);
		// Performs PCA
		myPCA.run();

		
		// System.out.println("EigenValues: " +myPCA.getEigenValues());
		
		DarkFieldEllipsoid myEllipsoid = new DarkFieldEllipsoid(scatterMatrix, scatterCoef, myPCA.getEigenValues(), myPCA.getEigenVectors());
		
		SimpleVector constrainedCoef = myEllipsoid.calculateSquaredProjectedCoefficients();
		
		return constrainedCoef;
		
	}
	
	
	/**
	 * Calculates the scatter direction of every voxel in the volume
	 * We use the ClassDarkFieldPCA to calculate the directions
	 * @return
	 */
	public DarkFieldTensorClass calcFiberOrientations( ){
		
		tensorClass = new DarkFieldTensorClass(imgSizeX, imgSizeY, imgSizeZ, darkFieldVolume.getSpacing(),darkFieldVolume.getOrigin());
		
		for(int x = 0; x < imgSizeX; x++){
			for(int y = 0; y <imgSizeY; y++){
				for(int z = 0; z <imgSizeZ; z++){
					// Calculate the scatter direction at given voxel
					DarkFieldPCA myPCA = getPCA_Result(x, y, z);
					// Add fiber direction to the ScatterDirectionObject
					tensorClass.setData(x,y,z,myPCA);
				} // End loop z
			} // End loop y
		} // End loop z
		
		return tensorClass;
		
	}
	
	/**
	 * Calculates the scatter direction at voxel element (x,y,z)
	 * @param x
	 * @param y
	 * @param z
	 * @return Scatter Direction
	 */
	private DarkFieldPCA getPCA_Result(int x, int y, int z){
		
		// Initializes the PCA objects
		SimpleVector scatterCoef = darkFieldVolume.getSimpleVectorAtIndex(x, y, z);
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterMatrix,scatterCoef);
		// Performs PCA
		myPCA.run();
		// Extract Scatter Direction as smallest component of the ellipsoid.
	
		return myPCA;
	}

	

	
	
}
