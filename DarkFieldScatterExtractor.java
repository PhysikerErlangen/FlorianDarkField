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
	 * Calculates the scatter direction of every voxel in the volume
	 * We use the ClassDarkFieldPCA to calculate the directions
	 * @return
	 */
	public DarkFieldFiberDirectionClass calcScatterDirections( ){
		
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
	public SimpleVector calcScatterDirection(int x, int y, int z){
		
		// Initializes the PCA objects
		SimpleVector scatterWeights = darkFieldVolume.getSimpleVectorAtIndex(x, y, z);
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterMatrix,scatterWeights);
		// Performs PCA
		myPCA.run();
		// Extract Scatter Direction as smallest component of the ellipsoid.
		
		double th = 1E-10;
		
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
