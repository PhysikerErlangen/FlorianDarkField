package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import edu.stanford.rsl.conrad.data.numeric.Grid4D;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldFiberDirectionClass extends Grid4D{

	public DarkFieldFiberDirectionClass(int imgSizeX, int imgSizeY,int imgSizeZ){
		super(imgSizeX, imgSizeY, imgSizeZ, 3);
	}
	
	@Override
	public void show(){
		
	}
	
	public void setFiberDirection(int x, int y, int z, SimpleVector direction){
		
		for(int i = 0; i < 3; i++){

			super.setAtIndex(x, y, z, i, (float) direction.getElement(i));
			
		}
	}
	
}
