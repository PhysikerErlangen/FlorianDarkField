// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 10st, 2014
//
package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.ImagePlus;
import edu.stanford.rsl.science.darkfield.FlorianDarkField.DarkField3DSinogram;


public class ImageToSinogram3D {
	
	
//	public static Grid2D imagePlusToGrid(DarkField3DSinogram img,int k){
//		Grid2D grid = new Grid2D(img.getSize()[1], img.getSize()[0]/2+1);
//		for(int j=0;j<img.getSize()[1];j++)
//			for(int i=0;i<img.getSize()[0]/2+1;i++){
//
//				grid.setAtIndex(j, i, img.getAtIndex(i, j, k));
//			}
//		grid.setSpacing(img.getSpacing()[1], img.getSpacing()[0]);
//
//
//		return grid;
//	}


	public static DarkField3DSinogram imagePlusToImagePlus3D(ImagePlus img){
		int maxThetaIndex =  img.getNSlices();
		int maxU = img.getWidth();
		int maxV = img.getHeight();
		
		DarkField3DSinogram grid = new DarkField3DSinogram(maxU,maxV,maxThetaIndex);

		for(int theta = 1; theta < maxThetaIndex +1; theta++){ // Start with Stack 1 so Matlab convention

			// An Image is actually saved as ONE vector, not a 2 dimensional image
			float[] sliceValues = ( float[]) img.getStack().getPixels(theta);
			
				for(int curU = 0; curU < maxU; curU++){
					for(int curV = 0; curV < maxV; curV++){
						
						int ind = curV*maxU+curU;
						
						// Do a interpolation of neighbor points of current value is 0
						if(sliceValues[ind]==0.f){
							// Add both neighbors up and take mean
							 if(ind == 0){
							sliceValues[ind]=sliceValues[ind+1];
							 } else if(ind == maxU*maxV -1){
							sliceValues[ind]=sliceValues[ind-1];
							 } else{
							sliceValues[ind]=(sliceValues[ind+1]+sliceValues[ind-1])/2.f;
							 } 
							 
							 if (sliceValues[ind] == 0.f){
								 sliceValues[ind] = 0.0001f; 
							 }
							 
						}
						
					double threshold = -Math.log( sliceValues[ind]); 
						
					if(  threshold < 0.f || Double.isNaN(threshold)){
						grid.setAtIndex(curU,curV,theta-1,0);
					}else{
						grid.setAtIndex(curU,curV,theta-1,(float) threshold);
						
				}
			}
			}
		}
		grid.setSpacing(1.d, 1.d,Math.PI*2/img.getNSlices());
		return grid;
	
	}
	
	// We need a different "conversion" algorithm for the absorption image
	// as DarkField Data has to through a logarithm first, while this is not valid for the
	// absorption image
	
	public static DarkField3DSinogram imagePlusToImagePlus3D_for_Absorption(ImagePlus img){
		
		int maxThetaIndex =  img.getNSlices();
		int maxU = img.getWidth();
		int maxV = img.getHeight();
		
		DarkField3DSinogram grid = new DarkField3DSinogram(maxU,maxV,maxThetaIndex);

		for(int theta = 1; theta < maxThetaIndex +1; theta++){ // Start with Stack 1 so Matlab convention

			// An Image is actually saved as ONE vector, not a 2 dimensional image
			float[] sliceValues = ( float[]) img.getStack().getPixels(theta);
			
				for(int curU = 0; curU < maxU; curU++){
					for(int curV = 0; curV < maxV; curV++){
						
						int ind = curV*maxU+curU;
						
						// Do a interpolation of neighbor points of current value is 0
						if(sliceValues[ind]==0.f){
							if(ind == 0){ // Check if one is at the left border
							sliceValues[ind]=sliceValues[ind+1];
							 } else if(ind == maxU*maxV -1){  // Check if one is at the right border
							sliceValues[ind]=sliceValues[ind-1];
							 } else{ // Add both neighbors up and take mean
							sliceValues[ind]=(sliceValues[ind+1]+sliceValues[ind-1])/2.f;
							 } 
							 if (sliceValues[ind] == 0.f){
								 sliceValues[ind] = 0.0001f; 
							 }
						}

					if(  sliceValues[ind]  < 0.f || Float.isNaN(sliceValues[ind]) ){
						grid.setAtIndex(curU,curV,theta-1,0);
					}else{
						grid.setAtIndex(curU,curV,theta-1,sliceValues[ind]);
					}
				}
			}
		}
		grid.setSpacing(1.d, 1.d,Math.PI*2/img.getNSlices());

		
		return grid;
	}
	}
//	public static Grid2D imagePlusToGridAMP(Grid3D img,int k){
//		Grid2D grid = new Grid2D(img.getSize()[1], img.getSize()[0]);
//		for(int j=0;j<img.getSize()[1];j++)
//			for(int i=0;i<img.getSize()[0];i++){
//				if(img.getAtIndex(i, j, k)>0.014)
//					grid.setAtIndex(j, i, img.getAtIndex(i, j, k));
//				else
//					grid.setAtIndex(j, i, 0);
//			}
//		grid.setSpacing(img.getSpacing()[1], img.getSpacing()[0]);
//
//
//		return grid;
//	}

	
	






