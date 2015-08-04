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
		
		fiberDirections = new DarkFieldFiberDirectionClass(imgSizeX, imgSizeY, imgSizeZ);
		
	}
	
	public void calcScatterDirections( ){
		
		for(int x = 0; x < imgSizeX; x++){
			for(int y = 0; y <imgSizeY; y++){
				for(int z = 0; z <imgSizeZ; z++){
					// Initializes the PCA objects
					SimpleVector scatterWeights = darkFieldVolume.getSimpleVectorAtIndex(x, y, z);
					DarkFieldPCA myPCA = new DarkFieldPCA(scatterMatrix,scatterWeights);
					// Performs PCA
					myPCA.run();
					// Extract Scatter Direction as smallest component of the ellipsoid.
					SimpleVector fiberDir = myPCA.getEigenVectors().getCol(2);
					
					fiberDirections.setFiberDirection(x, y, z, fiberDir);
				} // End loop z
			} // End loop y
		} // End loop z
		
		
		
	}
	
	
}
