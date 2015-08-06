package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;

import java.io.File;

/**
 * @author schiffers
 * DarkFieldReconPipeLine is the whole DarkField Reconstruction PipeLine
 * 
 * This Method basically contains 2 major functions:
 * 1. reconstructMaskForZeroConstraint
 * 	Deals with the reconstruction of the absorption data
 *  for calculating a binary mask, which can be used for a zero constraint
 *  during the DarkField reconstruction process
 * 
 * 2. reconstructDarkFieldVolume
 *  Is dealing with the DarkField Reconstruction
 *  A prior calculated reconstruction mask can be used. 
 *
 */
public class DarkFieldReconPipeline {

	
	Configuration config1;
	Configuration config2;
	
	private File fileDCI1;
	private File fileDCI2;
	private File fileAMP1;
	private File fileAMP2; //TODO Is not used at the current state
	
	boolean debug = true;

	
	/**
	 *  Reconstructed dark field volume
	 */
	private DarkField3DTensorVolume reconDarkField;
	
	/**
	 * @return the reconDarkField
	 */
	public DarkField3DTensorVolume getReconDarkField() {
		return reconDarkField;
	}


	/**
	 * Mask used for zero constrained enforcement during
	 * dark field reconstruction.
	 * Can be set to null if no mask is desired
	 */
	private DarkField3DTensorVolume reconMask;
	
	/**
	 * @return the reconMask
	 */
	public DarkField3DTensorVolume getReconMask() {
		return reconMask;
	}

	/**
	 * @param reconMask the reconMask to set
	 */
	public void setReconMask(DarkField3DTensorVolume reconMask) {
		this.reconMask = reconMask;
	}
	
	/**
	 * contains the different scatter directions used for reconstruction
	 */
	private SimpleMatrix scatterDirMatrix;
	
	

	/**
	 * CONSTRUCTOR
	 */
	
	/**
	 * @param fileNameConfig1
	 * @param fileNameDCI1
	 * @param fileNameDCI2
	 * @param fileNameAMP1
	 * @param fileNameAMP2
	 */
	public DarkFieldReconPipeline(String fileNameConfig1, String fileNameDCI1, String fileNameDCI2, 
			String fileNameAMP1, String fileNameAMP2){
		Configuration config1 = Configuration.loadConfiguration(fileNameConfig1);
		Configuration config2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		config2.getGeometry().setRotationAxis(rotationAxis2);
		initData(config1, config2, fileNameDCI1, fileNameDCI2, fileNameAMP1, fileNameAMP2);
		
	}
	
	/**
	 * @param fileNameConfig1
	 * @param fileNameConfig2
	 * @param fileNameDCI1
	 * @param fileNameDCI2
	 * @param fileNameAMP1
	 * @param fileNameAMP2
	 */
	public DarkFieldReconPipeline(String fileNameConfig1, String fileNameConfig2, String fileNameDCI1, String fileNameDCI2, 
			String fileNameAMP1, String fileNameAMP2){
		Configuration config1 = Configuration.loadConfiguration(fileNameConfig1);
		Configuration config2 = Configuration.loadConfiguration(fileNameConfig2);
		initData(config1, config2, fileNameDCI1, fileNameDCI2, fileNameAMP1, fileNameAMP2);
	}

	
	/**
	 * INITILIAZATION OF DATA NEEDED FOR THE RECONSTRUCTION OF THE DARKFIELD
	 * @param config1
	 * @param config2
	 * @param fileNameDCI1
	 * @param fileNameDCI2
	 * @param fileNameAMP1
	 * @param fileNameAMP2
	 */
	private void initData(Configuration config1, Configuration config2, 
			String fileNameDCI1, String fileNameDCI2, 
			String fileNameAMP1, String fileNameAMP2){
		
		new ImageJ();
		
		this.config1 = config1;
		this.config2 = config2;
		
		if((fileNameAMP1==null)||(fileNameAMP2==null)){
		fileAMP1 = null;			
		fileAMP2 = null;
			
		} else{
			this.fileAMP1 = new File(fileNameAMP1);
		this.fileAMP2 = new File(fileNameAMP2);
		}
		
		
		this.fileDCI1 = new File(fileNameDCI1);
		this.fileDCI2 = new File(fileNameDCI2);
		
		// Initialize the reconstruction mask to be null
		reconMask = null;
		
	}
	
	/**
	 ******************************************************************************
	 * METHODS
	 ****************************************************************************** 
	 */
	

	/**
	 * OVERLOADED FUNCTION THAT NEVERS SAVES THE DATA
	 * @param th_lower
	 * @param th_higher
	 */
	public void reconstructMaskForZeroConstraint(float th_lower, float th_higher){
		reconstructMaskForZeroConstraint(th_lower, th_higher, false, false);
	}

   /**
    * Calculates the 3D-Binary Mask which is used for zero Constraint 
	* in Dark-Field Tensor Reconstruction.
	* 
	*  The thresholds th_lower and th_higher define the threshold for 
	*  creating the mask out of the absorption volume. 
    * @param th_lower - lower threshold for binary mask
    * @param th_higher - upper threshold for binary mask
    * @param saveAMP - Boolean if AMP Image should be saved
    * @param saveMask - Boolean if Mask should be saved
    */
public void reconstructMaskForZeroConstraint(float th_lower, float th_higher, boolean saveAMP, boolean saveMask){
	   
	   /**
		 * INITILIAZATION OF ABSORPTION DATA
		 */
	   
		// Load absorption image of orientation 1
		ImagePlus imgAMP1 = IJ.openImage(fileAMP1.getAbsolutePath());
		DarkField3DSinogram sinoAMP1   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP1);

		/*
		 * JUST USE ONE ABSORPTION VOLUME FOR THE MOMENT
		 */
		
//		// Load absorption image of orientation 2
//		ImagePlus imgAMP2 = IJ.openImage(fileNameAMP2);
//		DarkField3DSinogram sinoAMP2   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP2);		
		
		if(debug) System.out.println("Load Parallel Beam Reconstruction Pipeline.");
		// Initialize the Parallel Beam Absorption Reconstruction
		DarkFieldAbsorptionRecon3D parallellBeamRecon = new DarkFieldAbsorptionRecon3D(config1);
		
		// Reconstruct the Absorption volume later used for zero constraint
				parallellBeamRecon.reconstructAbsorptionVolume(sinoAMP1);
		
		if(debug){
			System.out.println("Absorption Volume was reconstructed.");
		}
		
		// Save absorption reconstruction into a tif file
		if(saveAMP){
		   String filePath = fileAMP1.getParent() + "\\AMP_recon_volume.tif";
		   String volumeName = "Reconstruction of Absorption volume";
		   parallellBeamRecon.getReconAMP().write3DTensorToImage(filePath, volumeName);
		}
		
		// Create Mask out of absorption reconstruction
		reconMask = parallellBeamRecon.createMaskByBinaryThresholding(th_lower, th_higher);
		
		// Save absoprtion mask into a tif file 
		if(saveMask){
			String filePath = fileAMP1.getParent() + "\\AMP_mask_volume.tif";
			String volumeName = "Masked used for zero constraint";
			parallellBeamRecon.getMyMask().write3DTensorToImage(filePath, volumeName);
		}
		
		if(debug){
			System.out.println("Absorption mask was created and saved.");
		}
		
		// Dereference non used Data manually
		parallellBeamRecon = null;
		imgAMP1 = null;
		// imgAMP2 = null;
   }
   
   
	/**
	 * @param numScatterVectors
	 * @param maxIt
	 * @param stepSize
	 */
	public void reconstructDarkFieldVolume(int numScatterVectors, int maxIt, float stepSize, boolean saveDarkField){
		
		/*
		 * INITILIAZATION OF ABSORPTION DATA
		 */
	   
		// Load dark field image of orientation 1
		ImagePlus imgDCI1 = IJ.openImage(fileDCI1.getPath());
		// Readout DarkFieldSinogram out of ImagePlus Image
		DarkField3DSinogram sinoDCI1   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI1);
		
		// Load dark field image of orientation 2
		ImagePlus imgDCI2 = IJ.openImage(fileDCI2.getPath());
		// Readout DarkFieldSinogram out of ImagePlus Image
		DarkField3DSinogram sinoDCI2   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI2);		
		
		// Dereference unused data automatically
		imgDCI1 = null;
		imgDCI2 = null;
		
		// Initialize the GradientSolver3D
		// Can even deal with reconMask == null, as GradientSolver knows how to deal with this.
		// If mask should be used, it should be created prior to execution of this method
		GradientSolverTensor3D gradientSolver = new GradientSolverTensor3D(config1, config2, sinoDCI1, sinoDCI2, stepSize, maxIt, numScatterVectors,reconMask,reconMask);
		
		// Save the scatter directions into a matrix
		scatterDirMatrix =  gradientSolver.getScatterDirectionsAsMatrix();
		
		// Execute the gradient decent
		reconDarkField = gradientSolver.Gradient3D();

		if(saveDarkField){
			String filePath = fileDCI1.getParent() + "\\DCI_volume.tif";
			String volumeName = "Reconstructed DarkField Volume";
			reconDarkField.write3DTensorToImage(filePath, volumeName);
		}
	}
	

	/**
	 * @param saveVTK - if true saves the fiber Direcitons into a vtk file.
	 */
	public void calculateFiberOrientations(boolean saveVTK){
	
		calculateFiberOrientations(saveVTK, reconDarkField, scatterDirMatrix, fileDCI1);
	
	}
	
	

	/**
	 * @param saveVTK
	 * @param reconDarkField
	 * @param scatterDirMatrix
	 * @param fiberDirectionClass
	 * @param fileDCI1
	 */
	public static void calculateFiberOrientations(boolean saveVTK, DarkField3DTensorVolume reconDarkField, SimpleMatrix scatterDirMatrix, File fileDCI1){
		
		assert(reconDarkField != null) : new Exception("Reconstructed is NULL and needs to be calculated first.");
		assert(scatterDirMatrix != null) : new Exception("Scatter Dir Matrix needs to be calculated first.");
		
		DarkFieldScatterExtractor scatterDirExtractor = new DarkFieldScatterExtractor(reconDarkField, scatterDirMatrix);
		
		DarkFieldFiberDirectionClass fiberDirectionClass = new DarkFieldFiberDirectionClass
				(reconDarkField.imgSizeX, reconDarkField.imgSizeY, reconDarkField.imgSizeZ,
						reconDarkField.getSpacing(), reconDarkField.getOrigin());
		
		fiberDirectionClass = scatterDirExtractor.calcScatterDirections();
		
		
		if(saveVTK){
			String pathFiberVTK = fileDCI1.getParent() + "\\fiberDirections.vtk";
			fiberDirectionClass.writeToVectorField(pathFiberVTK);
		}
	}
	
	
}
