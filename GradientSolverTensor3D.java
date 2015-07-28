// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;

// Specialized backprojector and projector methods are required for solving the system with gradient decent
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ParallelDarkFieldBackprojector3DTensor; 
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ParallelDarkFieldProjector3DTensor;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;
import edu.stanford.rsl.conrad.utils.Configuration;


public class GradientSolverTensor3D extends DarkFieldTensorGeometry{
	
	
	 Configuration configuration1;
	 Configuration configuration2;
			
	DarkField3DSinogram darkFieldSinogram1;
	DarkField3DSinogram darkFieldSinogram2;
	
	float stepSize;
	
	int maxIt;
		
	DarkFieldScatterCoef scatterCoef1;
	DarkFieldScatterCoef scatterCoef2;
	
	
	// Constructor of GradientSolver3D
	public GradientSolverTensor3D( Configuration configuration1, Configuration configuration2, DarkField3DSinogram darkFieldSinogram1, 	DarkField3DSinogram darkFieldSinogram2, 
			float stepSize, int maxIt, int numScatterVectors){

		
		// Open super operator of geometry class
		super(configuration1,numScatterVectors);
		
		this.numScatterVectors = numScatterVectors;
		
		this.stepSize = stepSize;
		this.maxIt = maxIt;
		
		this.darkFieldSinogram1 = darkFieldSinogram1;
		this.darkFieldSinogram2 = darkFieldSinogram2;
	
		this.configuration1 = configuration1;
		this.configuration2 = configuration2;
		
	
		// Create instances of both scatter coef classes
		// One for each direction
		scatterCoef1 = new DarkFieldScatterCoef(configuration1,numScatterVectors);
		scatterCoef2 = new DarkFieldScatterCoef(configuration2,numScatterVectors);
	}
	
	/*
	 * Gradient 3D implements the gradient decent algorithm described in book of "zeng"  
	 */
	
	public DarkField3DTensorVolume Gradient3D() throws Exception{
		
		boolean debug = true;
		
		boolean reconVertical = true;
		boolean reconHorizontal = false;
		
		// Initialize to be constructed volume
		DarkField3DTensorVolume reconImage = new DarkField3DTensorVolume(
				imgSizeX,imgSizeY,imgSizeZ,numScatterVectors,getSpacing(),getOrigin());

		
		// Create instance of the backprojector
		ParallelDarkFieldBackprojector3DTensor backProjector1 = new ParallelDarkFieldBackprojector3DTensor(configuration1, scatterCoef1);
		ParallelDarkFieldBackprojector3DTensor backProjector2 = new ParallelDarkFieldBackprojector3DTensor(configuration2, scatterCoef2);
		
		// Create instance of the projector
		ParallelDarkFieldProjector3DTensor projector1 = new ParallelDarkFieldProjector3DTensor(configuration1, scatterCoef1);
		ParallelDarkFieldProjector3DTensor projector2 = new ParallelDarkFieldProjector3DTensor(configuration2,scatterCoef2);

		System.out.println("--------------------------------------------------------"  + maxIt);
		System.out.println("Start Gradient Decent Algorithm. Number of iteration is: " + maxIt);
		System.out.println("--------------------------------------------------------"  + maxIt);
		
		System.out.println("Stepsize: " +stepSize);
		
		
		// Iterate over all iterations
		for( int it = 0; it < maxIt; it++){
			
			long startTime = System.currentTimeMillis();
			
			System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
			System.out.println("Current Iteration: " + it);
			
			/* Gradient Decent
			1. Calculate difference between Measurement and Reconstruction by projection
			2. Backproject this difference onto the reconstruction (volume)
			3. Add this backprojected volume ontop of your current reconstruction
			*/
			
			DarkField3DSinogram projectionSinogram1 = null;
			DarkField3DSinogram projectionSinogram2 = null;
			
			if(it == 0){
				projectionSinogram1 = new DarkField3DSinogram(maxU_index, maxV_index, maxTheta_index);
				projectionSinogram2 = new DarkField3DSinogram(maxU_index, maxV_index, maxTheta_index);
			
			} else{

				// Calculate projection of first reconstruction
				if(reconVertical){
				projectionSinogram1 = projector1.projectPixelDriven(reconImage);
				}
				if(reconHorizontal){
				projectionSinogram2 = projector2.projectPixelDriven(reconImage);
				}
			
				if(debug){
					//projectionSinogram1.showSinogram("Sinogram 1 at it: " +it);
					//projectionSinogram2.showSinogram("Sinogram 2 at it: " +it);
					System.out.println("End projection of current reconstruction (It: " + it +")");
					System.out.println("Start calculation of difference between observation and reconstruction (It: " + it +")");
				}
				
			}
	
			DarkField3DSinogram differenceSinogram1 = null,differenceSinogram2 = null;
			
			// Calculate difference between observation and current projection
			if(reconVertical){
			differenceSinogram1 = DarkField3DSinogram.sub(projectionSinogram1,darkFieldSinogram1);
			differenceSinogram1.showSinogram("Sinogram1 of differences at iteration: " +it);

			// Backprojection difference between observation (Sinogram) and current iteration
			DarkField3DTensorVolume backProjectionDifference1 = backProjector1.backprojectPixelDriven(differenceSinogram1);
			// First multiply this with the gradient step size
			backProjectionDifference1.multiply(stepSize);
			reconImage.sub(backProjectionDifference1);
			}
			
			
			if(reconHorizontal){
			differenceSinogram2 = DarkField3DSinogram.sub(projectionSinogram2,darkFieldSinogram2);
			differenceSinogram2.showSinogram("Sinogram2 of differences at iteration: " +it);
			// Backprojection difference between observation (Sinogram) and current iteration
			DarkField3DTensorVolume backProjectionDifference2 = backProjector2.backprojectPixelDriven(differenceSinogram2);

			// Apply gradient step by adding difference on top of current reconstruction
			backProjectionDifference2.multiply(stepSize);
			reconImage.sub(backProjectionDifference2);
			}
			
//			if(debug){
//				backProjectionDifference1.show("Backprojection1 of difference at it:" +it);
//				backProjectionDifference2.show("Backprojection2 of difference at it:" +it);
//			}
//			
			

			reconImage.show("reconstruction at iteration: " +it);
			
			long endTime = System.currentTimeMillis();
			long deltaT = endTime - startTime;
			System.out.println("Gradient step completed in " + deltaT + "ms");
			
			}
		
		// Return the reconstruction result
		return reconImage;
		
	}	
}
