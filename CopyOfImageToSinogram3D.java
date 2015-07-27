// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers July 10st, 2014
//
package edu.stanford.rsl.science.darkfield.Florian;

import ij.ImagePlus;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.science.darkfield.Florian.DarkField3DSinogram;


public class CopyOfImageToSinogram3D {
	
	
	public static Grid2D imagePlusToGrid(DarkField3DSinogram img,int k){
		Grid2D grid = new Grid2D(img.getSize()[1], img.getSize()[0]/2+1);
		for(int j=0;j<img.getSize()[1];j++)
			for(int i=0;i<img.getSize()[0]/2+1;i++){

				grid.setAtIndex(j, i, img.getAtIndex(i, j, k));
			}
		grid.setSpacing(img.getSpacing()[1], img.getSpacing()[0]);


		return grid;
	}


	public static DarkField3DSinogram imagePlusToImagePlus3D(ImagePlus img){
		DarkField3DSinogram grid = new DarkField3DSinogram(img.getNSlices(), img.getWidth(), img.getHeight());
		for(int i=1;i<img.getNSlices()+1;i++)
			for(int j=0;j<img.getWidth();j++)
				for(int k=0;k<img.getHeight();k++){
					float[]floats=(float[])img.getStack().getPixels(i);
					if(floats[k*img.getWidth()+j]==0.f)
						floats[k*img.getWidth()+j]=(floats[k*img.getWidth()+j+1]+floats[k*img.getWidth()+j-1])/2.f;
					if(-Math.log( floats[k*img.getWidth()+j]) < 0.f)
						grid.setAtIndex(i-1,j, k,0);
					else
						grid.setAtIndex(i-1, j,k,(float) -Math.log( floats[k*img.getWidth()+j]) );
				}


		grid.setSpacing(Math.PI*2/img.getNSlices(),1.d, 1.d);


		return grid;
	}
	public static Grid2D imagePlusToGridAMP(DarkField3DSinogram img,int k){
		Grid2D grid = new Grid2D(img.getSize()[1], img.getSize()[0]);
		for(int j=0;j<img.getSize()[1];j++)
			for(int i=0;i<img.getSize()[0];i++){
				if(img.getAtIndex(i, j, k)>0.014)
					grid.setAtIndex(j, i, img.getAtIndex(i, j, k));
				else
					grid.setAtIndex(j, i, 0);
			}
		grid.setSpacing(img.getSpacing()[1], img.getSpacing()[0]);


		return grid;
	}

	public static DarkField3DSinogram imagePlusToImagePlus3DAMP(ImagePlus img){
		DarkField3DSinogram grid = new DarkField3DSinogram(img.getNSlices(), img.getWidth(), img.getHeight());
		for(int i=1;i<img.getNSlices()+1;i++)
			for(int j=0;j<img.getWidth();j++)
				for(int k=0;k<img.getHeight();k++){
					float[]floats=(float[])img.getStack().getPixels(i);
					grid.setAtIndex(i-1, j,k, floats[k*img.getWidth()+j] );
				}


		grid.setSpacing(Math.PI*2/img.getNSlices(),1.d,1.d);


		return grid;
	}
}




