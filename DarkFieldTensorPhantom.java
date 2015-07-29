package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.util.ArrayList;

import weka.core.pmml.jaxbbindings.GaussianDistribution;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import edu.stanford.rsl.apps.gui.opengl.PointCloudViewer;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.conrad.utils.FileUtil;
import edu.stanford.rsl.conrad.utils.RegKeys;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.DarkFieldPhantom;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.ImageToGrid3D;
import edu.stanford.rsl.science.darkfield.iterative.GradientSolver;
import edu.stanford.rsl.science.darkfield.reconstruction.WriteToVectorField;
import edu.stanford.rsl.tutorial.basics.MHDImageLoader;
import edu.stanford.rsl.tutorial.basics.PointCloudMaker;



public class DarkFieldTensorPhantom  extends  DarkFieldTensorGeometry  {

	
	ArrayList<DarkField3DSinogram> sinogramList;
	ArrayList<ImagePlus> projectionList;
	
	DarkField3DTensorVolume phantom;
	
	public static void main (String [] args) throws Exception{
	
		String fileNameConfig1 = "C:\\Users\\schiffers\\workspace\\Configurations\\PhantomConfigurationFlorian-001-axis.xml";
		// Load configuration wooden case
		Configuration config = Configuration.loadConfiguration(fileNameConfig1);
		Configuration Configuration2 = Configuration.loadConfiguration(fileNameConfig1);
		// Reset rotation axis for Config2
		SimpleVector rotationAxis2 = new SimpleVector(0.0d,1.0d,0.0);
		Configuration2.getGeometry().setRotationAxis(rotationAxis2);
		
		
		DarkFieldTensorPhantom myPhantom = new DarkFieldTensorPhantom(config);

		new ImageJ();
		
		myPhantom.phantom.show();
		
		myPhantom.calculateDarkFieldProjection(config,Configuration2);
		
		myPhantom.getDarkFieldSinogram(0).show();
		myPhantom.getDarkFieldSinogram(1).show();
		
		myPhantom.getDarkFieldSinogram(0).showSinogram();
		myPhantom.getDarkFieldSinogram(1).showSinogram();
		
		
		System.out.println("Phaton created");
	}



	
	public DarkFieldTensorPhantom(Configuration config){
		this(config,3);
	}
	

	
	public DarkFieldTensorPhantom(Configuration config, int numChannels){

		// Call super constructor of TensorGeometry
		super(config,numChannels);
		
		phantom = new DarkField3DTensorVolume(this.imgSizeX,this.imgSizeY, this.imgSizeZ, numChannels,this.getSpacing(),this.getOrigin());
		
		sinogramList = new ArrayList<DarkField3DSinogram>();
		
		calcPhantom();
	}
	
	
	public void calcPhantom(){
		
		
		// USE THREE DIFFERENT SCATTER DIRECTION
		
		SimpleVector scatterDir1 = new SimpleVector(1,0,0);
		SimpleVector scatterDir2 = new SimpleVector(0,1,0);
		SimpleVector scatterDir3 = new SimpleVector(0,0,1);
		
		
		int a = (int)( 0.2*imgSizeX);
		int b = (int) (0.8*imgSizeY);
		
		int a2 = (int)( 0.4*imgSizeX);
		int b2 = (int) (0.6*imgSizeY);
		
		int a3 = (int)( 0.3*imgSizeX);
		int b3 = (int) (0.7*imgSizeY);
		
		System.out.println("a:" + a + "b: " + b);
		
		for(int x = a; x < b; x ++){
			for(int y = a2; y < b2; y ++){
				for(int z = a3; z < b3; z ++){

					// Creates a phantom of 3 Layers, each Layer contains voxel with a constant scatter direction!
					
					phantom.setDarkFieldScatterTensor(x, y, z, scatterDir1);
					
//					if((double)z <= a ){
//							phantom.setDarkFieldScatterTensor(x, y, z, scatterDir1);
//					}else if( (double)z < b){
//							phantom.setDarkFieldScatterTensor(x, y, z, scatterDir2);
//					}else{
//							phantom.setDarkFieldScatterTensor(x, y, z, scatterDir3);
//					}

				} // END X 
		} // END Y
	} // END Z
	

} // End calcPhantom()
	
	
	
	
	public DarkField3DSinogram getDarkFieldSinogram(int index){
		return sinogramList.get(index);
	}
	
	public void calculateDarkFieldProjection(Configuration config1, Configuration config2){
		
		DarkFieldScatterCoef scatterCoef1 = new DarkFieldScatterCoef(config1,phantom.getNumberOfChannels());
		DarkFieldScatterCoef scatterCoef2 = new DarkFieldScatterCoef(config2,phantom.getNumberOfChannels());
		
		ParallelDarkFieldProjector3DTensor projector1 = new
				ParallelDarkFieldProjector3DTensor(config1, scatterCoef1);
		
		ParallelDarkFieldProjector3DTensor projector2 = new
				ParallelDarkFieldProjector3DTensor(config2, scatterCoef2);
		
		//DarkField3DSinogram sino1 = projector1.projectRayDriven(phantom);
		//DarkField3DSinogram sino2 = projector2.projectRayDriven(phantom);
		
		DarkField3DSinogram sino1 = projector1.projectPixelDriven(phantom);
		DarkField3DSinogram sino2 = projector2.projectPixelDriven(phantom);
		
		
		
		sinogramList.add(sino1);
		sinogramList.add(sino2);
		
	}
	
}