// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 16st, 2015
//

package edu.stanford.rsl.science.darkfield.Florian;


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

	public void setMaxThetaIndex(int maxThetaIndex) {
		this.maxThetaIndex = maxThetaIndex;
	}

	public DarkField3DSinogram(int maxU, int maxV,int maxThetaIndex){
		super(maxU, maxV, maxThetaIndex,true);
		
		this.maxU = maxU;
		this.maxV = maxV;
		this.maxThetaIndex = maxThetaIndex;
		
	}
	

	// Multiplies the whole grid with a given factor
	public void multiply(float factor){
		
		for(int u = 0; u <this.getSize()[0]; u++){
			for(int v = 0; v <this.getSize()[1]; v++){
				for(int theta = 0; theta <this.getSize()[0]; theta++){
					multiplyAtIndex(u,v,theta,factor);
				} // End loop theta
			} // End loop v
		} // End loop u
	}
	
	
	public void show(String title){
		ImagePlus img = ImageUtil.wrapGrid3D(this, title);
		 
		img.show();

	}
	
	public void show(){
		show("");
	}
	
	public void showSinogram(){
	showSinogram("");
	}
	
	public void showSinogram(String title){
		
		Grid3D mySinoGrid = new Grid3D(maxThetaIndex,maxU,maxV);
		
		for(int curTheta = 0; curTheta < maxThetaIndex; curTheta++){ // Start with Stack 1 so Matlab convention
			for(int curU = 0; curU < maxU; curU++){
				for(int curV = 0; curV < maxV; curV++){
					mySinoGrid.setAtIndex(curTheta,curU,curV,getAtIndex(curU, curV, curTheta));
				}
			}
		}

		ImageUtil.wrapGrid3D(mySinoGrid, title).show();
		
	
	}
	
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
	
	public boolean checkForNaN(int curU, int curV, int curTheta){
		float myVal = getAtIndex(curU,curV,curTheta);
		if (Float.isNaN(myVal)){
			return true;
		}
		return false;
	}
	
	
	public void sub( DarkField3DSinogram B) throws Exception {

		// Check for inconstiency (different dimensions)
		if(getSize()[0]!=B.getSize()[0]&&getSize()[1]!=B.getSize()[1]&&getSize()[2]!=B.getSize()[2]){
		System.out.println("Dimensions do not match2");
		throw new Exception("Dimension do not match!");
		}
		
		for(int u = 0; u <getMaxU(); u++){
			for(int v = 0; v <getMaxV(); v++){
				for(int theta = 0; theta <getMaxThetaIndex(); theta++){
				this.addAtIndex(u, v, theta, - B.getAtIndex(u, v, theta));
				} // End loop theta
			} // End loop v
		} // End loop u
		
	}
	

	
	public static DarkField3DSinogram sub(DarkField3DSinogram A, DarkField3DSinogram B) throws Exception{

		// Check for inconsistency (different dimensions)
		if(A.getSize()[0]!=B.getSize()[0]&&A.getSize()[1]!=B.getSize()[1]){
		System.out.println("Dimensions do not match in substraction of 2 sinograms");
		throw new Exception("Dimension do not match in substraction of 2 sinograms!");
		}
		
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
	

}
