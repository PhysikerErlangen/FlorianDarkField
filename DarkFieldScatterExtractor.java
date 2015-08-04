package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldScatterExtractor {

	
	DarkField3DTensorVolume darkFieldVolume;
	SimpleMatrix scatterMatrix;
	
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
	
	public static void calcScatterDirections( ){
		
		for(int x = 0; x <myVolume.imgSizeX; x++){
			for(int y = 0; y <myVolume.imgSizeY; y++){
				for(int z = 0; z <myVolume.imgSizeZ; z++){
					// Initializes the PCA objects
					DarkFieldPCA myPCA = new DarkFieldPCA(scatterMatrix,scatterWeights);
					// Performs PCA
					myPCA.run();
					// Extract Scatter Direction as smallest component of the ellipsoid.
					
				} // End loop z
			} // End loop y
		} // End loop z
		
		
		
	}
	
	
}
