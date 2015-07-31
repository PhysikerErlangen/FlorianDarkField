//This code was developed in a collaboration with ECAP, Erlangen, Germany.
//This part of the code is not to be published under GPL before Oct 31st 2017.
//author@ Florian Schiffers July 1st, 2015
//
package edu.stanford.rsl.science.darkfield.FlorianDarkField;
import ij.IJ;
import ij.ImagePlus;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.data.numeric.Grid4D;
import edu.stanford.rsl.conrad.utils.ImageUtil;


public class DarkField3DTensorVolume extends DarkFieldGrid3DTensor{

	@SuppressWarnings("unused")
	private String title;
	

	public DarkField3DTensorVolume(int imageSizeX,int imageSizeY, int imageSizeZ, int numChannels,double[] spacing, double[] origin){
		this(imageSizeX, imageSizeY, imageSizeZ, numChannels, spacing, origin, "");
	}

	public DarkField3DTensorVolume(int imageSizeX,int imageSizeY, int imageSizeZ, int numChannels,double[] spacing, double[] origin, String title) {
		
		// Call superconstructor of DarkFieldGrid3DTensor
		super(imageSizeX, imageSizeY, imageSizeZ,numChannels);
		
		// Set spacing of box
		setSpacing(spacing);
		// Set origin of 3D Image Box
		setOrigin(origin);
		// set Title
		setTitle(title);
	}

	
	@Override
	public Grid3D getChannel(int c){
		Grid3D myGrid =  getMultichannelData().getSubGrid(c);
		myGrid.setSpacing(getSpacing());
		myGrid.setOrigin(getOrigin());
		return myGrid;
	}
	
	public void setTitle(String t) {
		this.title = t;
	}
	
 // masks a volume with a given mask. Careful, has to have same dimension
	public void maskWithVolume(DarkField3DTensorVolume mask){
		
		// If there's no mask just return and do nothing
		if(mask == null)
			return;
	
		int[] mySize = getSize();
		
		for(int x = 0; x < mySize[0]; x++){
			for(int y = 0; y < mySize[1]; y++){
				for(int z = 0; z < mySize[2]; z++){
					// The mask may only have ! one ! channel as it is as a mask
					if(mask.getPixelValue(x, y, z, 0) == 0f){
					this.setAtIndex(x, y, z, 0f);
					}
			}
		}
	}
	}		
	
	
	
	public void sub( DarkField3DTensorVolume B) throws Exception{

		// Check for inconstiency (different dimensions)
		if(getSize()[0]!=B.getSize()[0]&&getSize()[1]!=B.getSize()[1]&&getSize()[2]!=B.getSize()[2]){
		System.out.println("Dimensions do not match2");

		System.out.println(getSize());
		System.out.println(B.getSize());
		
		throw new Exception("Dimension do not match in substraction of 2 tensor volumes!");
		
		}
		
		for(int x = 0; x <this.getSize()[0]; x++){
			for(int y = 0; y <this.getSize()[1]; y++){
				for(int z = 0; z <this.getSize()[0]; z++){
					
				this.subAtDarkFieldScatterTensor(x, y, z, B.getVectorAtIndex(x, y, z));
				
				} // End loop z
			} // End loop y
		} // End loop z
		
	}
	
	
	
	
//	
//	@Override
//	public void show(){
//		show("");
//		
//	}
//	
//	@Override
//	public void show(String reconTitle){
//		
//		this.getMultichannelData().show();
//		
//		
//		for(int channel = 0; channel < this.getNumberOfChannels(); channel++){
//			
//			Grid3D curGrid = this.getChannel(channel);
//			
//			String title = reconTitle + "Channel Nr. "+channel;
//			ImagePlus myImg = ImageUtil.wrapGrid3D(curGrid, title);
//			
//			myImg.show();
//		
//			IJ.run(myImg, "Volume Viewer", "");
//			// curGrid.show();
//		}
//			
//			
//	}
	
}
