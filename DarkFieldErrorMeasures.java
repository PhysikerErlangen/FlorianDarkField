package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.data.numeric.NumericGridOperator;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldErrorMeasures {

	
	public static enum DarkFieldNormType{
		NORM_L1,
		NORM_L2,
		NORM_MAX
	}
	
	
	
	/**
	 * Calculates the  angular distance of all reconstructed fiber directions
	 * @param A
	 * @param B
	 * @return
	 */
	public static Grid3D errorAngularDistanceGrid(DarkFieldTensorClass A, DarkFieldTensorClass B){
		
		// Check for inconsistency (different dimensions)
				assert(A.imgSizeX == B.imgSizeX
						&&A.imgSizeY == B.imgSizeY
						&&A.imgSizeZ == B.imgSizeZ
						): new Exception("Dimension of data is wrong.");
					
	Grid3D angularDiffGrid = new Grid3D(A.imgSizeX, A.imgSizeY, A.imgSizeZ);
	
	DarkFieldVectorField dirA = A.getFiberDirection();
	DarkFieldVectorField dirB = B.getFiberDirection();
	
	 for(int x = 0; x < A.imgSizeX; x ++){
		 for(int y = 0; y < A.imgSizeY; y++){
			 for(int z = 0; z < A.imgSizeZ; z++){
				 
				 SimpleVector vecA  = dirA.getSimpleVectorAtIndex(x, y, z).normalizedL2();
				 SimpleVector vecB  = dirB.getSimpleVectorAtIndex(x, y, z).normalizedL2();
				 
				 // Calculate inner product
				 double inner = SimpleOperators.multiplyInnerProd(vecA,vecB);
				 
				 // Take absolute value
				 inner = Math.abs(inner);
				 
				 /*
				  * Calculate angle and set it to grid
				  */
				 angularDiffGrid.setAtIndex(x,y,z,(float)Math.acos(inner));
				 
			 }
		 }
	 }
		return angularDiffGrid;
	}
	
	
	/**
	 * Calculates the average angular distance of all reconstructed fiber directions
	 * @param A
	 * @param B
	 * @return
	 */
	public static double errorAngularDistance(DarkFieldTensorClass A, DarkFieldTensorClass B){
		
		Grid3D diff = errorAngularDistanceGrid(A,B);
		double norm = NumericGridOperator.getInstance().normL1(diff);
		return norm;
	}
	
	
	/**
	 * Normalized residual norms 
	 * (see. Vogel - Constrained ... eq. 21 on page 13)
	 * @param A - Normalized in respect to A
	 * @param B
	 * @param normType
	 * @return
	 */
	public static double errorSinogam(DarkField3DSinogram A, DarkField3DSinogram B,DarkFieldNormType normType){
		DarkField3DSinogram diff = 	DarkField3DSinogram.sub(A,B);
		if(normType == DarkFieldNormType.NORM_L1){
		double residualNorm = diff.norm1();
		double normalisedResidualNorm = residualNorm/A.norm1();
		return normalisedResidualNorm;
		} else if(normType == DarkFieldNormType.NORM_L2){
		double residualNorm = diff.norm2();
		double normalisedResidualNorm = residualNorm/A.norm2();
		return normalisedResidualNorm;
		}else{
			return Double.NaN;
		}
	}
	

	/**
	 * @param diff
	 * @param normFactor
	 * @param normType
	 * @return
	 */
	public static double errorSinogam(DarkField3DSinogram diff, double normFactor,DarkFieldNormType normType){
		if(normType == DarkFieldNormType.NORM_L1){
		double residualNorm = diff.norm1();
		double normalisedResidualNorm = residualNorm/normFactor;
		return normalisedResidualNorm;
		} else if(normType == DarkFieldNormType.NORM_L2){
		double residualNorm = diff.norm2();
		double normalisedResidualNorm = residualNorm/normFactor;
		return normalisedResidualNorm;
		}else{
			return Double.NaN;
		}
	}
	

	/**
	 * Normalised mean updates
	 * (see eq. 22 of Vogel - constrained XTT)
	 * @param A
	 * @param B - Normalized in respect to B
	 * @return
	 */
	public static double errorDarkFieldCoefficients(DarkField3DTensorVolume A, DarkField3DTensorVolume B,DarkFieldNormType normType ){
		DarkField3DTensorVolume diff = DarkField3DTensorVolume.sub(A, B);
		
		if(normType == DarkFieldNormType.NORM_L1){
		return diff.norm1()/B.norm1();
		} else if(normType == DarkFieldNormType.NORM_L2){
		return diff.norm2()/B.norm2();
		}else{
			return Double.NaN;
		}
		
	}
	
	
	/**
	 * @param folder
	 * @param fileName
	 * @param data
	 */
	public static void writeErrorToTxt(File folder, String fileName, SimpleMatrix data){
		
		  try {
	            
	            File newTextFile = new File(folder.getParent() + "\\" + fileName);
	            FileWriter fw = new FileWriter(newTextFile);
	            
	            for(int row = 0; row < data.getRows(); row++){
	            	
	            	SimpleVector rowVec = data.getRow(row);
	            	
	            	for( int col = 0; col < data.getCols(); col++){
	            		fw.write(rowVec.getElement(col) + " ");	   
	            	}
	            	fw.write(System.getProperty( "line.separator" ));
	            }
	            fw.close();
	        } catch (IOException iox) {
	            //do stuff with exception
	            iox.printStackTrace();
	        }
		
	}
	 
	
}
