// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;

// Specialized backprojector and projector methods are required for solving the system with gradient decent
import java.io.File;

import edu.stanford.rsl.science.darkfield.FlorianDarkField.ParallelDarkFieldBackprojector3DTensor;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ParallelDarkFieldProjector3DTensor;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;
import edu.stanford.rsl.science.darkfield.iterative.OpMath;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.data.numeric.iterators.NumericPointwiseIteratorND;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.utils.Configuration;

public class GradientSolverTensor3D extends DarkFieldTensorGeometry {

	boolean debug = false;

	boolean reconVertical = true;
	boolean reconHorizontal = false;

	boolean writeVtkInEveryStep = true;
	
	// private Configuration configuration1;
	// private Configuration configuration2;

	private DarkField3DSinogram darkFieldSinogram1;
	private DarkField3DSinogram darkFieldSinogram2;

	private float stepSize;

	private int maxIt;

	DarkFieldScatterCoef scatterCoef1;
	DarkFieldScatterCoef scatterCoef2;

	// MASKING Images for zero constraint
	DarkField3DTensorVolume reconAMP1;
	DarkField3DTensorVolume reconAMP2;

	private ParallelDarkFieldBackprojector3DTensor backProjector1;
	private ParallelDarkFieldBackprojector3DTensor backProjector2;

	private ParallelDarkFieldProjector3DTensor projector1;
	private ParallelDarkFieldProjector3DTensor projector2;

	DarkField3DTensorVolume reconImage;

	private File pathToSaveVtk;
	
	
	/**
	 * @param configuration1
	 * @param configuration2
	 * @param darkFieldSinogram1
	 * @param darkFieldSinogram2
	 * @param stepSize
	 * @param maxIt
	 * @param numScatterVectors
	 * @param pathToSaveVtk
	 */
	public GradientSolverTensor3D(Configuration configuration1,
			Configuration configuration2,
			DarkField3DSinogram darkFieldSinogram1,
			DarkField3DSinogram darkFieldSinogram2, float stepSize, int maxIt,
			int numScatterVectors, File pathToSaveVtk) {
		this(configuration1, configuration2, darkFieldSinogram1,
				darkFieldSinogram2, stepSize, maxIt, numScatterVectors, pathToSaveVtk, null,
				null);
	}

	/**
	 * @param configuration1
	 * @param configuration2
	 * @param darkFieldSinogram1
	 * @param darkFieldSinogram2
	 * @param stepSize
	 * @param maxIt
	 * @param numScatterVectors
	 * @param pathToSaveVtk
	 * @param reconAMP1
	 * @param reconAMP2
	 */
	public GradientSolverTensor3D(Configuration configuration1,
			Configuration configuration2,
			DarkField3DSinogram darkFieldSinogram1,
			DarkField3DSinogram darkFieldSinogram2, float stepSize, int maxIt,
			int numScatterVectors,  File pathToSaveVtk, DarkField3DTensorVolume reconAMP1,
			DarkField3DTensorVolume reconAMP2) {

		// Open super operator of geometry class
		super(configuration1, numScatterVectors);

		this.numScatterVectors = numScatterVectors;
		
		this.stepSize = stepSize;
		this.maxIt = maxIt;

		this.darkFieldSinogram1 = darkFieldSinogram1;
		this.darkFieldSinogram2 = darkFieldSinogram2;

		// this.configuration1 = configuration1;
		// this.configuration2 = configuration2;

		this.reconAMP1 = reconAMP1;
		this.reconAMP2 = reconAMP2;
		
		this.pathToSaveVtk = pathToSaveVtk;

		// Create instances of both scatter coef classes
		// One for each direction
		scatterCoef1 = new DarkFieldScatterCoef(configuration1,
				numScatterVectors);
		scatterCoef2 = new DarkFieldScatterCoef(configuration2,
				numScatterVectors);

		// Create instance of the backprojector
		backProjector1 = new ParallelDarkFieldBackprojector3DTensor(
				configuration1, scatterCoef1);
		backProjector2 = new ParallelDarkFieldBackprojector3DTensor(
				configuration2, scatterCoef2);

		// Create instance of the projector
		projector1 = new ParallelDarkFieldProjector3DTensor(configuration1,
				scatterCoef1);
		projector2 = new ParallelDarkFieldProjector3DTensor(configuration2,
				scatterCoef2);

	}

	/**
	 * Gradient 3D implements the gradient decent algorithm described in book of
	 * "zeng"
	 * 
	 * @return reconstructed DarkField Tensor volume
	 */
	public DarkField3DTensorVolume Gradient3D() {

		debug = true;
		
		if(darkFieldSinogram1==null){
		reconVertical = false;
		}else{
		reconVertical = true;
		}
		
		if(darkFieldSinogram1==null){
		reconHorizontal = false;
		}else{
			reconHorizontal = true;
		}
		
		writeVtkInEveryStep = true;


		// Initialize to be constructed volume
		reconImage = new DarkField3DTensorVolume(imgSizeX, imgSizeY, imgSizeZ,
				numScatterVectors, getSpacing(), getOrigin());

		reconImage.show("Current Iteration of Reconstructed Volume");
		//reconImage.showComponents();

		System.out
				.println("--------------------------------------------------------"
						+ maxIt);
		System.out
				.println("Start Gradient Decent Algorithm. Number of iteration is: "
						+ maxIt);
		System.out
				.println("--------------------------------------------------------"
						+ maxIt);

		System.out.println("Stepsize: " + stepSize);

		/*
		 *  Iterate over all iterations
		 *  and perform on gradient step in each iteration
		 */
		
		for (int it = 0; it < maxIt; it++) {

			long startTime = System.currentTimeMillis();

			System.out
					.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
			System.out.println("Current Iteration: " + it);

			doGradientStep(it);

			if(writeVtkInEveryStep){
				
				String myName = "ReconstructedVolumeIter_"+it;
				reconImage.saveFiberOrientations(pathToSaveVtk, myName);
				
			}
			
			long endTime = System.currentTimeMillis();
			long deltaT = endTime - startTime;
			System.out.println("Gradient step completed in " + deltaT
					+ "ms, It: " + it);

		}

		// Return the reconstruction result
		return reconImage;

	}

	/**
	 * Perform one step of the gradient
	 */
	private void doGradientStep(int it) {

		long startTime = System.currentTimeMillis();

		System.out
				.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		System.out.println("Current Iteration: " + it);

		long endTime = System.currentTimeMillis();
		long deltaT = endTime - startTime;
		System.out.println("Gradient step completed in " + deltaT + "ms, It: "
				+ it);

		/*
		 * Gradient Decent 1. Calculate difference between Measurement and
		 * Reconstruction by projection 2. Backproject this difference onto the
		 * reconstruction (volume) 3. Add this backprojected volume ontop of
		 * your current reconstruction
		 */

		DarkField3DSinogram projectionSinogram1 = null;
		DarkField3DSinogram projectionSinogram2 = null;

		int numElements = maxU_index * maxV_index * maxTheta_index;
		double error1 = 0, error2 = 0;

		if (debug)
			System.out.println("Start projection of current Iteration.");

		if (it == 0) {
			
			if (reconVertical) {
				projectionSinogram1 = new DarkField3DSinogram(maxU_index,
						maxV_index, maxTheta_index);
			}if (reconHorizontal) {
				projectionSinogram2 = new DarkField3DSinogram(maxU_index,
						maxV_index, maxTheta_index);
			}
		} else {
			// Calculate projection of first reconstruction
			if (reconVertical) {
				projectionSinogram1 = projector1.projectPixelDriven(reconImage);
			} if (reconHorizontal) {
				projectionSinogram2 = projector2.projectPixelDriven(reconImage);
			}

			if (debug) {
				// projectionSinogram1.showSinogram("Sinogram 1 at it: " +it);
				// projectionSinogram2.showSinogram("Sinogram 2 at it: " +it);
				System.out
						.println("End projection of current reconstruction (It: "
								+ it + ")");
				System.out
						.println("Start calculation of difference between observation and reconstruction (It: "
								+ it + ")");
			}

		}

		if (debug)
			System.out.println("End projection of current Iteration.");

		DarkField3DSinogram differenceSinogram1 = null, differenceSinogram2 = null;

		// Calculate difference between observation and current projection
		if (reconVertical) { // Reconstruct vertical trajectory

			if (debug)
				System.out.println("Start reconstruction of Trajectory 1.");

			differenceSinogram1 = DarkField3DSinogram.sub(projectionSinogram1,
					darkFieldSinogram1);
			error1 = (float) OpMath.norm2(differenceSinogram1) / numElements;
			// Backprojection difference between observation (Sinogram) and
			// current iteration

			if (debug)
				System.out
						.println("Start Backprojection of Differences of Trajector 1.");
			DarkField3DTensorVolume backProjectionDifference1 = backProjector1
					.backprojectPixelDriven(differenceSinogram1);
			if (debug)
				System.out
						.println("End Backprojection of Differences of Trajector 1.");
			// First multiply this with the gradient step size
			// differenceSinogram1.showSinogram("Difference of Sinograms");

			// backProjectionDifference1.

			backProjectionDifference1.multiply(stepSize);
			backProjectionDifference1.maskWithVolume(reconAMP1);

			reconImage.sub(backProjectionDifference1);

			if (debug)
				System.out.println("End reconstruction of Trajectory 1.");
		}

		if (reconHorizontal) {// Reconstruct horizontal trajectory

			if (debug)
				System.out.println("Start reconstruction of Trajectory 2.");

			differenceSinogram2 = DarkField3DSinogram.sub(projectionSinogram2,
					darkFieldSinogram2);
			error2 = (float) OpMath.norm2(differenceSinogram2) / numElements;
			// differenceSinogram2.showSinogram("Sinogram2 of differences at iteration: "
			// +it);
			if (debug)
				System.out
						.println("Start Backprojection of Differences of Trajector 2.");
			// Backprojection difference between observation (Sinogram) and
			// current iteration
			DarkField3DTensorVolume backProjectionDifference2 = backProjector2
					.backprojectPixelDriven(differenceSinogram2);
			// Apply gradient step by adding difference on top of current
			// reconstruction
			backProjectionDifference2.multiply(stepSize);
			backProjectionDifference2.maskWithVolume(reconAMP2);
			reconImage.sub(backProjectionDifference2);

			if (debug)
				System.out.println("End reconstruction of Trajectory 2.");
		}

		double totalError = error1 + error2;
		System.out.println("Error (Difference of Sinograms): " + totalError);

	}

}
