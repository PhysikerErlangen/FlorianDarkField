package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldScatterExtractor {

	
	private DarkField3DTensorVolume darkFieldVolume;
	private SimpleMatrix scatterMatrix;
	
	private DarkFieldFiberDirectionClass fiberDirections;
	
	private int imgSizeX;
	private int imgSizeY;
	private int imgSizeZ;
	
	
	/**
	 * Threshold that checks, if 3 component of eigenvalues is too small
	 * If 3 component is too small, don't consider it as a fiber orientation
	 * and ignore it
	 */
	double th = 1E-10;
	
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
	public DarkFieldFiberDirectionClass calcFiberOrientations( ){
		
		fiberDirections = new DarkFieldFiberDirectionClass(imgSizeX, imgSizeY, imgSizeZ, darkFieldVolume.getSpacing(),darkFieldVolume.getOrigin());
		
		for(int x = 0; x < imgSizeX; x++){
			for(int y = 0; y <imgSizeY; y++){
				for(int z = 0; z <imgSizeZ; z++){
					// Calculate the scatter direction at given voxel
					SimpleVector fiberDir = calcScatterDirection(x, y, z);
					// Add fiber direction to the ScatterDirectionObject
					fiberDirections.setFiberDirection(x, y, z, fiberDir);
				} // End loop z
			} // End loop y
		} // End loop z
		
		return fiberDirections;
		
	}
	
	/**
	 * Calculates the scatter direction at voxel element (x,y,z)
	 * @param x
	 * @param y
	 * @param z
	 * @return Scatter Direction
	 */
	private SimpleVector calcScatterDirection(int x, int y, int z){
		
		// Initializes the PCA objects
		SimpleVector scatterCoef = darkFieldVolume.getSimpleVectorAtIndex(x, y, z);
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterMatrix,scatterCoef);
		// Performs PCA
		myPCA.run();
		// Extract Scatter Direction as smallest component of the ellipsoid.
		
		
		
		SimpleVector fiberDir;
		if(myPCA.getEigenValues().getElement(2)<th){
			fiberDir = new SimpleVector(3);
		}else{
			fiberDir = myPCA.getEigenVectors().getCol(2).normalizedL2();
			fiberDir.multiplyBy(myPCA.getEigenValues().getElement(2));
		}
		
		
		return fiberDir;
	}

	
}
