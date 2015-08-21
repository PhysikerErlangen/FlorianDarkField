// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 16st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;


import java.io.File;

import ij.IJ;
import ij.ImagePlus;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.utils.ImageUtil;


public class DarkField3DSinogram extends Grid3D {


	private int maxU;
	private int maxV;
	private int maxThetaIndex;
	
 // DarkField3DSinogramm works with float values
		
	public int getMaxU() {
		return maxU;
	}

	public void setMaxU(int maxU) {
		this.maxU = maxU;
	}

	public int getMaxV() {
		return maxV;
	}

	public void setMaxV(int maxV) {
		this.maxV = maxV;
	}

	public int getMaxThetaIndex() {
		return maxThetaIndex;
	}

	/**
	 * @param maxThetaIndex
	 */
	public void setMaxThetaIndex(int maxThetaIndex) {
		this.maxThetaIndex = maxThetaIndex;
	}

	/**
	 * @param maxU
	 * @param maxV
	 * @param maxThetaIndex
	 */
	public DarkField3DSinogram(int maxU, int maxV,int maxThetaIndex){
		super(maxU, maxV, maxThetaIndex,true);
		
		this.maxU = maxU;
		this.maxV = maxV;
		this.maxThetaIndex = maxThetaIndex;
		
	}
	@Override
	public float getAtIndex(int curU, int curV, int thetaIdx){
		return super.getAtIndex(curU, curV, thetaIdx);
	}
	
	/**
	 * Multiplies the whole grid with a given factor
	 * @param factor
	 */
	public void multiply(float factor){
		
		for(int u = 0; u <this.getSize()[0]; u++){
			for(int v = 0; v <this.getSize()[1]; v++){
				for(int theta = 0; theta <this.getSize()[0]; theta++){
					multiplyAtIndex(u,v,theta,factor);
					
				} // End loop theta
			} // End loop v
		} // End loop u
	}
	
	
	@Override
	public void show(String title){
		ImagePlus img = ImageUtil.wrapGrid3D(this, title);
				
		img.show();
		
		IJ.run(img, "Enhance Contrast", "saturated=0.35");

	}
	
	@Override
	public void show(){
		show("");
	}
	
	/**
	 * Displays sinogam with void title ""
	 */
	public void showSinogram(){
	showSinogram("");
	}
	
	/**
	 * Displays the sinogram.
	 * Has to create a new Grid which is why this might take some time.
	 * @param title is the title of the sinogram
	 */
	public void showSinogram(String title){

		ImagePlus img = ImageUtil.wrapGrid3D(this, title);
	
		IJ.run(img, "Reslice [/]...", "output=1.000 start=Top rotate avoid");
		
		
//		Grid3D mySinoGrid = new Grid3D(maxThetaIndex,maxU,maxV);
//		
//		double[] mySpacing = this.getSpacing();
//		mySinoGrid.setSpacing(mySpacing[2],mySpacing[1],mySpacing[0]);
//		
//		double[] myOrigin = this.getOrigin();
//		
//		mySinoGrid.setOrigin(myOrigin[2],myOrigin[1],myOrigin[0]);
//		
//		for(int curTheta = 0; curTheta < maxThetaIndex; curTheta++){ // Start with Stack 1 so Matlab convention
//			for(int curU = 0; curU < maxU; curU++){
//				for(int curV = 0; curV < maxV; curV++){
//					mySinoGrid.setAtIndex(curTheta,curU,curV,getAtIndex(curU, curV, curTheta));
//				}
//			}
//		}
//
//		ImagePlus img = ImageUtil.wrapGrid3D(mySinoGrid, title);
//		
//		img.show();
	}
	
	/**
	 * Checks if one of all Elements in the sinogram is NaN
	 * @return
	 */
	public boolean checkForNan(){
		
		for(int curTheta = 0; curTheta < maxThetaIndex; curTheta++){ // Start with Stack 1 so Matlab convention
			for(int curU = 0; curU < maxU; curU++){
				for(int curV = 0; curV < maxV; curV++){
					float myVal = getAtIndex(curU, curV, curTheta);
					// System.out.println(myVal);
					if (Float.isNaN(myVal)){
						System.out.println("First NaN at curU = " + curU + " curV = " +curV + " theta " +curTheta);
						return true;
					}
				}
			}
		}
		
		return false;
		
	}
	
	/**
	 * Checks if value of sinogram is naN
	 * @param curU
	 * @param curV
	 * @param curTheta
	 * @return true if is NaN, false if it is not a NaN
	 */
	public boolean checkForNaN(int curU, int curV, int curTheta){
		float myVal = getAtIndex(curU,curV,curTheta);
		if (Float.isNaN(myVal)){
			return true;
		}
		return false;
	}

	/**
	 * @return
	 */
	public double norm1(){
		return norm1(this);
	}
	
	
	/**
	 * @param sinogram
	 * @return
	 */
	public static double norm1(DarkField3DSinogram sinogram) {
		double result=0.d;
		
		for(int thetaIdx=0;thetaIdx<sinogram.maxThetaIndex;thetaIdx++){
			for(int curU=0;curU<sinogram.maxU;curU++){
				for(int curV=0;curV<sinogram.maxV; curV++){
					result+= Math.abs(sinogram.getAtIndex(curU, curV, thetaIdx));
				}
			}
		}
		return result;
	}

	
	
	/**
	 * Subtracts B from current Sinogram
	 * @param B
	 */
	public void sub( DarkField3DSinogram B) {

		// Check for inconsistency (different dimensions)
				assert( getSize()[0]==B.getSize()[0]
						&&getSize()[1]==B.getSize()[1]
						&&getSize()[2]!=B.getSize()[2])
						: new Exception("Dimension of data is wrong.");
		
		for(int u = 0; u <getMaxU(); u++){
			for(int v = 0; v <getMaxV(); v++){
				for(int theta = 0; theta <getMaxThetaIndex(); theta++){
				this.addAtIndex(u, v, theta, - B.getAtIndex(u, v, theta));
				} // End loop theta
			} // End loop v
		} // End loop u
		
	}
	

	
	/**
	 * @param A - Sinogram 1
	 * @param B - Sinogram 2
	 * @return Sinogram subtracted A - B
	 */
	public static DarkField3DSinogram sub(DarkField3DSinogram A, DarkField3DSinogram B){

		// Check for inconsistency (different dimensions)
		assert( (A.getSize()[0]==B.getSize()[0])
				&&
				A.getSize()[1]==B.getSize()[1]
						)
				: new Exception("Dimension of data is wrong.");
		
		// Create new instance of a DarkField3DSinogram with same dimensions
		DarkField3DSinogram C= new DarkField3DSinogram(A.getSize()[0],A.getSize()[1],A.getSize()[2]);
		
		for(int u = 0; u <A.getMaxU(); u++){
			for(int v = 0; v <A.getMaxV(); v++){
				for(int theta = 0; theta <A.getMaxThetaIndex(); theta++){
				C.setAtIndex(u, v, theta, A.getAtIndex(u, v, theta)- B.getAtIndex(u, v, theta));
				} // End loop theta
			} // End loop v
		} // End loop u
		return C;
	}
	
	
	/**
	 * writes the file to a given file specified by filePath
	 * The name of the Image is not further defined, so default ""
	 * @param folderPath - Path where Volume should be saved
	 */
	public void writeDarkFieldSinogramToImage(String folderPath){
		writeDarkFieldSinogramToImage(folderPath, "DarkFieldSinogram");
	}
	
	
	/**
	 * @param folderPath - "C\\User\\MyFolder\\"
	 * @param sinogramName - "MySinogramTitle"
	 * FileName is created out of sinogram Name
	 */
	public void writeDarkFieldSinogramToImage(String folderPath, String sinogramName){
		String fileName = sinogramName + ".vtk";
		writeDarkFieldSinogramToImage(new File(folderPath),fileName,sinogramName);
		}
	
	
	/**
	 * @param folderPath   - e.g. "C\\User\\MyFolder\\"
	 * @param fileName     - e.g. "DarkFieldSinogram1.tif"
	 * @param sinogramName - e.g. "MyDarkFieldSignal1"
	 */
	public void writeDarkFieldSinogramToImage(File folderPath, String fileName,  String sinogramName){
		ImagePlus img = ImageUtil.wrapGrid3D(this, sinogramName);
		String totalPath = folderPath.getParent() + "\\" + fileName;
		IJ.save(img,totalPath);
	}
	
	
	
}
