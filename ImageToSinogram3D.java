// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 10st, 2014
//
package edu.stanford.rsl.science.darkfield.Florian;

import ij.ImagePlus;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.science.darkfield.Florian.DarkField3DSinogram;


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
						
						// Do a interpolation of neighbor points of current value is 0
						if(sliceValues[curV*maxU+curU]==0.f){
							// Add both neighbors up and take mean
							sliceValues[curV*maxU+curU]=(sliceValues[curV*maxV+curU+1]+sliceValues[curV*maxU+curU-1])/2.f;
						}
					if( -Math.log( sliceValues[curV*maxU+curU]) < 0.f){
						grid.setAtIndex(curU,curV,theta-1,0);
					}else{
						grid.setAtIndex(curU,curV,theta-1,(float) -Math.log( sliceValues[curV*maxU+curU]));
						
				}
			}
			}
		}
		grid.setSpacing(1.d, 1.d,Math.PI*2/img.getNSlices());
		return grid;
	
	}
	
//	
//	public static Grid2D imagePlusToGridAMP(DarkField3DSinogram img,int k){
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

	

}




