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
	
	private File saveFolder;
	
	boolean debug = true;

	
	/**
	 *  Reconstructed dark field volume
	 */
	private DarkField3DTensorVolume reconDarkField;
	
	/**
	 * @param reconDarkField the reconDarkField to set
	 */
	public void setReconDarkField(DarkField3DTensorVolume reconDarkField) {
		this.reconDarkField = reconDarkField;
	}

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
			String fileNameSaveFolder){
		Configuration config1 = Configuration.loadConfiguration(fileNameConfig1);
		Configuration config2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		config2.getGeometry().setRotationAxis(rotationAxis2);
		initData(config1, config2, fileNameDCI1, fileNameDCI2,fileNameSaveFolder);
		
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
			String fileNameSaveFolder){
		Configuration config1 = Configuration.loadConfiguration(fileNameConfig1);
		Configuration config2 = Configuration.loadConfiguration(fileNameConfig2);
		initData(config1, config2, fileNameDCI1, fileNameDCI2,fileNameSaveFolder);
	}

	
	
	/**
	 * @param fileNameConfig1
	 * @param fileNameConfig2
	 * @param fileNameDCI1
	 * @param fileNameDCI2
	 * @param fileNameAMP1
	 * @param fileNameAMP2
	 */
	public DarkFieldReconPipeline(Configuration config1, Configuration config2,String fileNameSaveFolder){
		
		initData(config1, config2, null, null,fileNameSaveFolder);
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
			String fileNameSaveFolder){
		
		if(fileNameSaveFolder==null){
			saveFolder = null;	
		} else{
			saveFolder = new File(fileNameSaveFolder);
		}
		
		
		
		this.config1 = config1;
		this.config2 = config2;

		
		if(fileNameDCI1==null){
			fileDCI1 = null;			
			} 
		else{
				fileDCI1 = new File(fileNameDCI1);
			}
		
		if(fileNameDCI2==null){
			fileDCI2 = null;
				
			} else{
				fileDCI2 = new File(fileNameDCI2);
			}
		
	
		
		// Initialize the reconstruction mask to be null
		reconMask = null;
		
	}
	
	/**
	 ******************************************************************************
	 * METHODS
	 ****************************************************************************** 
	 */
	
	
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
	public void reconstructMaskForZeroConstraint(float th_lower, float th_higher, boolean saveAMP, boolean saveMask, String pathToAMP){
		   /**
			 * INITILIAZATION OF ABSORPTION DATA
			 */
		   
			// Load absorption image of orientation 1
			ImagePlus imgAMP1 = IJ.openImage(pathToAMP );
			DarkField3DSinogram sinoAMP1   = ImageToSinogram3D.imagePlusToImagePlus3D_for_Absorption(imgAMP1);

			reconstructMaskForZeroConstraint(th_lower, th_higher, saveAMP, saveMask,sinoAMP1);
	}
	
	

	/**
	 * OVERLOADED FUNCTION THAT NEVERS SAVES THE DATA
	 * @param th_lower
	 * @param th_higher
	 */
	public void reconstructMaskForZeroConstraint(float th_lower, float th_higher,String fileAMP1){
		reconstructMaskForZeroConstraint(th_lower, th_higher, false, false,fileAMP1);
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
	
public void reconstructMaskForZeroConstraint(float th_lower, float th_higher, boolean saveAMP, boolean saveMask, DarkField3DSinogram sinoAMP){


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
				parallellBeamRecon.reconstructAbsorptionVolume(sinoAMP);
		
		if(debug){
			System.out.println("Absorption Volume was reconstructed.");
			parallellBeamRecon.getReconAMP().show("Reconstruction of Absorption CT");
		}
		
		// Save absorption reconstruction into a tif file
		if(saveAMP){
		   String filePath = fileDCI1.getParent() + "\\AMP_recon_volume.tif";
		   String volumeName = "Reconstruction of Absorption volume";
		   parallellBeamRecon.getReconAMP().write3DTensorToImage(filePath, volumeName);
		}
		
		// Create Mask out of absorption reconstruction
		reconMask = parallellBeamRecon.createMaskByBinaryThresholding(th_lower, th_higher);
		
		// Save absoprtion mask into a tif file 
		if(saveMask){
			String filePath = fileDCI1.getParent() + "\\AMP_mask_volume.tif";
			String volumeName = "Masked used for zero constraint";
			parallellBeamRecon.getMyMask().write3DTensorToImage(filePath, volumeName);
		}
		
		if(debug){
			System.out.println("Absorption mask was created and saved.");
		}
		
   }
   

/**
 * Sinograms and read directly out of the file paths that are given by the constructor
 * These sinograms objects are then passed to the overloaded function of reconstructedDarkFieldVolume
 * This is done, so that you could also work with Phantom data, where only the sinogram is passed
 * and not the path to the .tif image itself
 * @param numScatterVectors
 * @param maxIt
 * @param stepSize
 * @param pathDarkField
 */
public void reconstructDarkFieldVolume(int numScatterVectors, int maxIt, float stepSize, File pathDarkField){
	
	
	/*
	 * INITILIAZATION OF ABSORPTION DATA
	 */
   
	
	DarkField3DSinogram sinoDCI1 = null;
	
	if(fileDCI1 != null){
		// Load dark field image of orientation 1
		ImagePlus imgDCI1 = IJ.openImage(fileDCI1.getPath());
		// Readout DarkFieldSinogram out of ImagePlus Image
		sinoDCI1   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI1);
		// Dereference unused data automatically
		imgDCI1 = null;
	}
	
	DarkField3DSinogram sinoDCI2 = null;
	
	if(fileDCI2 != null){
		// Load dark field image of orientation 1
		ImagePlus imgDCI2 = IJ.openImage(fileDCI2.getPath());
		// Readout DarkFieldSinogram out of ImagePlus Image
		sinoDCI2   = ImageToSinogram3D.imagePlusToImagePlus3D(imgDCI2);
		// Dereference unused data automatically
		imgDCI2 = null;
	}

	reconstructDarkFieldVolume(numScatterVectors, maxIt, stepSize, pathDarkField, sinoDCI1, sinoDCI2);
	
	
}


	/**
	 * @param numScatterVectors
	 * @param maxIt
	 * @param stepSize
	 */
	public void reconstructDarkFieldVolume(int numScatterVectors, int maxIt, float stepSize, File pathDarkField,DarkField3DSinogram sinoDCI1,DarkField3DSinogram sinoDCI2){
		
				
		// Initialize the GradientSolver3D
		// Can even deal with reconMask == null, as GradientSolver knows how to deal with this.
		// If mask should be used, it should be created prior to execution of this method
		GradientSolverTensor3D gradientSolver = new GradientSolverTensor3D(config1, config2, sinoDCI1, sinoDCI2, stepSize, maxIt, numScatterVectors,saveFolder, reconMask,reconMask);
		
		// Save the scatter directions into a matrix
		scatterDirMatrix =  DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);
		
		// Execute the gradient decent
		reconDarkField = gradientSolver.Gradient3D();

		if(pathDarkField!=null){
			String filePath = pathDarkField.getParent() + "\\DCI_volume.tif";
			String volumeName = "Reconstructed DarkField Volume";
			reconDarkField.write3DTensorToImage(filePath, volumeName);
		}
	}
	

	/**
	 * @param saveVTK - if true saves the fiber Direcitons into a vtk file.
	 */
	public DarkFieldFiberDirectionClass calculateFiberOrientations(boolean saveVTK){
	
		if(saveVTK == false)
			return calculateFiberOrientations(reconDarkField, scatterDirMatrix, null);
		else{
			return calculateFiberOrientations(reconDarkField, scatterDirMatrix, fileDCI1);	
		}
			
		
	
	}
	
	/**
	 * @param saveVTK - if true saves the fiber Direcitons into a vtk file.
	 */
	public DarkFieldFiberDirectionClass calculateFiberOrientations(File pathFile){
	
		return calculateFiberOrientations(reconDarkField, scatterDirMatrix, pathFile);
	
	}
	

	/**
	 * @param saveVTK
	 * @param reconDarkField
	 * @param scatterDirMatrix
	 * @param fiberDirectionClass
	 * @param pathFile
	 */
	public static DarkFieldFiberDirectionClass calculateFiberOrientations(DarkField3DTensorVolume reconDarkField, SimpleMatrix scatterDirMatrix, File pathFile){
		
		return calculateFiberOrientations(reconDarkField, scatterDirMatrix, pathFile,"fiberDirections");
		
		
	}
	


	/**
	 * @param reconDarkField
	 * @param scatterDirMatrix
	 * @return
	 */
	public static DarkFieldFiberDirectionClass calculateFiberOrientations(DarkField3DTensorVolume reconDarkField, SimpleMatrix scatterDirMatrix){
		
		return calculateFiberOrientations(reconDarkField, scatterDirMatrix, null,null);
		
		
	}
	
	/**
	 * @param reconDarkField
	 * @param scatterDirMatrix
	 * @param pathFile
	 * @param fileName - without file type, as .vtk is added automatically
	 */
	public static DarkFieldFiberDirectionClass calculateFiberOrientations(DarkField3DTensorVolume reconDarkField, SimpleMatrix scatterDirMatrix, File pathFile, String fileName){
		
		assert(reconDarkField != null) : new Exception("Reconstructed is NULL and needs to be calculated first.");
		assert(scatterDirMatrix != null) : new Exception("Scatter Dir Matrix needs to be calculated first.");
		
		DarkFieldScatterExtractor scatterDirExtractor = new DarkFieldScatterExtractor(reconDarkField, scatterDirMatrix);
		
		DarkFieldFiberDirectionClass fiberDirectionClass = new DarkFieldFiberDirectionClass
				(reconDarkField.imgSizeX, reconDarkField.imgSizeY, reconDarkField.imgSizeZ,
						reconDarkField.getSpacing(), reconDarkField.getOrigin());
		
		fiberDirectionClass = scatterDirExtractor.calcFiberOrientations();
		
		
		if(pathFile!=null){
			String pathFiberVTK = pathFile.getParent() + "\\ " + fileName + ".vtk";
			fiberDirectionClass.writeToVectorField(pathFiberVTK);
		}
		
		return fiberDirectionClass;
	}
	
	
}
