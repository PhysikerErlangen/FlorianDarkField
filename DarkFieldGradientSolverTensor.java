// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 1st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.gui.Plot;
import java.io.File;







import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkFieldErrorMeasures.DarkFieldNormType;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ParallelDarkFieldBackprojector3DTensor;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.ParallelDarkFieldProjector3DTensor;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DTensorVolume.TensorConstraintType;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;

public class DarkFieldGradientSolverTensor extends DarkFieldTensorGeometry {

	boolean debug = false;

	boolean reconVertical = true;
	boolean reconHorizontal = false;

	/**
	 * 
	 */
	boolean writeVtkInEveryStep = true;
	
	/**
	 * 
	 */
	private DarkField3DSinogram darkFieldSinogram1;
	private DarkField3DSinogram darkFieldSinogram2;

	/**
	 *  Stepsize for Gradient Descent
	 */
	private float stepSize;

	/**
	 * Maximum Number of Iterations
	 */
	private int maxIt;
	
	/**
	 * Contains error of reconstruction in all steps
	 */
	SimpleMatrix errorMat;

	/**
	 * Scatter weights used for Tensor reconstruction
	 */
	private DarkFieldScatterWeightsCalculator scatterWeights1;
	private DarkFieldScatterWeightsCalculator scatterWeights2;


	/**
	 * MASKING Images for zero constraint 
	 */
	private DarkField3DTensorVolume maskAMP1;
	private DarkField3DTensorVolume maskAMP2;

	private ParallelDarkFieldBackprojector3DTensor backProjector1;
	private ParallelDarkFieldBackprojector3DTensor backProjector2;

	private ParallelDarkFieldProjector3DTensor projector1;
	private ParallelDarkFieldProjector3DTensor projector2;

	/**
	 * Volume to reconstructed by Gradient Method
	 */
	private DarkField3DTensorVolume reconImage;

	/**
	 * @return the reconImage
	 */
	public DarkField3DTensorVolume getReconImage() {
		return reconImage;
	}

	/**
	 * @param reconImage the reconImage to set
	 */
	public void setReconImage(DarkField3DTensorVolume reconImage) {
		this.reconImage = reconImage;
	}





	private File pathToSaveVtk;
	
	private TensorConstraintType tensorConstraint = TensorConstraintType.HARD_CONSTRAINT;
	
	
	/**
	 * @return the tensorConstraint
	 */
	public TensorConstraintType getTensorConstraint() {
		return tensorConstraint;
	}

	/**
	 * @param tensorConstraint the tensorConstraint to set
	 */
	public void setTensorConstraint(TensorConstraintType tensorConstraint) {
		this.tensorConstraint = tensorConstraint;
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
	 */
	public DarkFieldGradientSolverTensor(Configuration configuration1,
			Configuration configuration2,
			DarkField3DSinogram darkFieldSinogram1,
			DarkField3DSinogram darkFieldSinogram2, float stepSize, int maxIt,
			int numScatterVectors, File pathToSaveVtk,TensorConstraintType type) {
		this(configuration1, configuration2, darkFieldSinogram1,
				darkFieldSinogram2, stepSize, maxIt, numScatterVectors, pathToSaveVtk, null,
				null,type);
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
	 * @param maskAMP1
	 * @param maskAMP2
	 */
	public DarkFieldGradientSolverTensor(Configuration configuration1,
			Configuration configuration2,
			DarkField3DSinogram darkFieldSinogram1,
			DarkField3DSinogram darkFieldSinogram2, float stepSize, int maxIt,
			int numScatterVectors,  File pathToSaveVtk, DarkField3DTensorVolume maskAMP1,
			DarkField3DTensorVolume maskAMP2,TensorConstraintType type) {

		// Open super operator of geometry class
		super(configuration1, numScatterVectors);

		this.numScatterVectors = numScatterVectors;
		
		this.stepSize = stepSize;
		this.maxIt = maxIt;

		this.darkFieldSinogram1 = darkFieldSinogram1;
		this.darkFieldSinogram2 = darkFieldSinogram2;

		// this.configuration1 = configuration1;
		// this.configuration2 = configuration2;

		this.maskAMP1 = maskAMP1;
		this.maskAMP2 = maskAMP2;
		
		this.pathToSaveVtk = pathToSaveVtk;
		
		/*
		 * If the darkFieldSinograms are null (i.e. they don't exist)
		 * do not reconstruct 
		 */
		if(darkFieldSinogram1==null){
		reconVertical = false;
		}else{
		reconVertical = true;
		}
		
		if(darkFieldSinogram2==null){
		reconHorizontal = false;
		}else{
			reconHorizontal = true;
		}
		
		

		/*
		 * Create instances of both scatter coef classes
		 * One for each direction
		 */
		scatterWeights1 = new DarkFieldScatterWeightsCalculator(configuration1,
				numScatterVectors);
		scatterWeights2 = new DarkFieldScatterWeightsCalculator(configuration2,
				numScatterVectors);

		/*
		 * Create instances of the BackProjector
		 */
		backProjector1 = new ParallelDarkFieldBackprojector3DTensor(
				configuration1, scatterWeights1);
		backProjector2 = new ParallelDarkFieldBackprojector3DTensor(
				configuration2, scatterWeights2);

		/*
		 * Create instances of the projectors
		 */
		projector1 = new ParallelDarkFieldProjector3DTensor(configuration1,
				scatterWeights1);
		projector2 = new ParallelDarkFieldProjector3DTensor(configuration2,
				scatterWeights2);
		
		setTensorConstraint(type);
		
		errorMat = new SimpleMatrix(maxIt,2);

	}

	/**
	 * Gradient 3D implements the gradient decent algorithm described in book of
	 * "zeng"
	 * 
	 * @return reconstructed DarkField Tensor volume
	 */
	public DarkField3DTensorVolume Gradient3D(boolean writeVtkInEveryStep) {

		debug = true;


		// Initialize to be constructed volume
		reconImage = new DarkField3DTensorVolume(imgSizeX, imgSizeY, imgSizeZ,
				numScatterVectors, getSpacing(), getOrigin());

		reconImage.show("Current Iteration of Reconstructed Volume");
		//reconImage.showComponents();

		System.out
				.println("--------------------------------------------------------"
						+ maxIt);
		System.out
				.println("Start Gradient Decent Algorithm. Number of maximum iteration is: "
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

			/**
			 * GRADIENT STEP FOR BOTH TRAJECTORIES
			 * AT ITERATION IT
			 */
			doGradientStep(it);

			
			/**
			 * CONSTRAINT ENFORCEMENT
			 */
			reconImage.enforceConstraint(tensorConstraint);
			
			
			
			if(writeVtkInEveryStep){
				System.out.println("Write Reconstruction Result at iteration = " +it);
				String myName = "ReconstructedVolumeIter_"+it;
				reconImage.saveFiberOrientations(pathToSaveVtk, myName);
				
			}
			
			long endTime = System.currentTimeMillis();
			long deltaT = endTime - startTime;
			System.out.println("Gradient step completed in " + deltaT
					+ "ms, It: " + it);

		}

		plotError(this.errorMat, this.maxIt);
		
		DarkFieldErrorMeasures.writeErrorToTxt(pathToSaveVtk, "error.txt", errorMat);
		
		// Return the reconstruction result
		return reconImage;

	}
		
	/**
	 * Perform one step of the gradient.
	 * Deals with both trajectories. 
	 * Reconstructed of each trajectory is handled in a submethod.
	 */
	private void doGradientStep(int it) {

		System.out
				.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		System.out.println("Current Iteration: " + it);

		/*
		 * Gradient Decent
		 * 1. Calculate difference between Measurement and reconstruction by projection
		 * 2. Backproject this difference onto the reconstruction (volume)
		 * 3. Add this backprojected volume on top of your current reconstruction
		 */

		double error = 0;
		
	
		if(reconVertical){
		double val =reconstructionTrajectory(TrajectoryType.VERTICAL,it);
		error = error +  val;
		}
		
		if(reconHorizontal){
		double val = reconstructionTrajectory(TrajectoryType.HORIZONTAL,it);
		error = error +  val;
		}	
	
		errorMat.setRowValue(it, new SimpleVector(it,error));
		
		plotError(errorMat, it+1);
		
		System.out.println("Error (Difference of Sinograms): " + error);

	}
	
	/**

	 */
	private void plotError(SimpleMatrix errorMatrix, int it){
		double[] data = errorMat.getCol(1).getSubVec(0, it ).copyAsDoubleArray();
		Plot myErrorPlot = VisualizationUtil.createPlot(data, "Error of Sinograms", "Iteration", "Error");
		myErrorPlot.show();
	}
	
	
	/**
	 * @param type
	 * @return
	 */
	private ParallelDarkFieldBackprojector3DTensor getBackProjector(TrajectoryType type){
	
		
		
		if(type == TrajectoryType.VERTICAL ){
			return backProjector1;
		} else if(type == TrajectoryType.HORIZONTAL ){
			return backProjector2;
		} else{
			return null;
		}
	}
	

	/**
	 * @param type
	 * @return
	 */
	private ParallelDarkFieldProjector3DTensor getProjector(TrajectoryType type){
		if(type == TrajectoryType.VERTICAL ){
			return projector1;
		} else if(type == TrajectoryType.HORIZONTAL ){
			return projector2;
		} else{
			return null;
		}
	}
	
	/**
	 * @param type
	 * @return
	 */
	private DarkField3DTensorVolume getMask(TrajectoryType type){
		if(type == TrajectoryType.VERTICAL ){
			return maskAMP1;
		} else if(type == TrajectoryType.HORIZONTAL ){
			return maskAMP2;
		} else{
			return null;
		}
	}
	
	/**
	 * @param type
	 * @return
	 */
	private DarkField3DSinogram getSinogram(TrajectoryType type){
		if(type == TrajectoryType.VERTICAL ){
			return darkFieldSinogram1;
		} else if(type == TrajectoryType.HORIZONTAL ){
			return darkFieldSinogram2;
		} else{
			return null;
		}
	}
	
	
	/**
	 * Reconstruction of one trajectory in Gradient Step
	 * @param darkFieldSinogram
	 * @param projector
	 * @param backprojector
	 * @param iteration
	 */
	private double reconstructionTrajectory(TrajectoryType type,
			int iteration){
		
		ParallelDarkFieldBackprojector3DTensor backProjector = getBackProjector(type);
		ParallelDarkFieldProjector3DTensor projector = getProjector(type);
		DarkField3DSinogram darkFieldSinogram = getSinogram(type);
		
		/*
		 * Calculate DarkFieldProjection of current reconstruction state
		 * First iteration handled differently, as projection signal (sinogram)
		 * of a 0-Volume is 0 anyway
		 */
		DarkField3DSinogram projectionSinogram;
		if(iteration == 0){
			projectionSinogram = new DarkField3DSinogram(maxU_index,
				maxV_index, maxTheta_index);
		} else{
			projectionSinogram = projector.projectPixelDriven(reconImage,getMask(type));	
		}
		
		/*
		 * Calculate difference between observed darkFieldSignal and reconstruction 
		 */
		DarkField3DSinogram differenceSinogram = 
				DarkField3DSinogram.sub(projectionSinogram,darkFieldSinogram);

		/*
		 * Backprojection of the projection difference between observation and reconstruction 
		 */
		DarkField3DTensorVolume diffVolume = backProjector
				.backprojectPixelDriven(differenceSinogram,getMask(type));

		/*
		 * Multiply diffVolume with stepSize of Gradient
		 */
		
		diffVolume.multiply(calcStepSize(stepSize, iteration));
		
		/*
		 * Mask the Image with the volume.
		 * MaskWithVolume handles null data
		 */
		//diffVolume.maskWithVolume(getMask(type));
		
		reconImage.sub(diffVolume);
		
		/*
		 * Calculate error between observed and reconstruction sinograms
		 */

		double error = DarkFieldErrorMeasures.errorSinogam(differenceSinogram,darkFieldSinogram.norm2(),DarkFieldNormType.NORM_L2);
		
		return error;
	}

	/**
	 * Calculates stepSize dependent on initial value and current iteration
	 * @param alpha
	 * @param it
	 */
	private float calcStepSize(float stepSize, int it) {
		return stepSize;
		// return  (float)( Math.pow(0.96, it)*stepSize);
	}
	
	
}
