package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid4D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.science.darkfield.darkfieldgrid.DarkFieldPhantom;

public class DarkFieldFiberDirectionClass extends Grid4D{

	// Defines dimension of volume box
	protected int imgSizeX;
	protected int imgSizeY;
	protected int imgSizeZ;
	
	

	/**
	 * @param imgSizeX
	 * @param imgSizeY
	 * @param imgSizeZ
	 * @param spacing_world
	 * @param origin_world
	 */
	public DarkFieldFiberDirectionClass(int imgSizeX, int imgSizeY,int imgSizeZ, double[] spacing_world, double[] origin_world){
		super(imgSizeX, imgSizeY, imgSizeZ, 3);
		
		this.imgSizeX = imgSizeX;
		this.imgSizeY = imgSizeY;
		this.imgSizeZ = imgSizeZ;
		
		// Set spacing of box
		setSpacing(spacing_world);
		// Set origin of 3D Image Box
		setOrigin(origin_world);
		
	}
	
	
	@Override
	public void show(){
		
	}
	
	/**
	 * @param x
	 * @param y
	 * @param z
	 * @param direction
	 */
	public void setFiberDirection(int x, int y, int z, SimpleVector direction){
		
		// Loop through every coordinate
		for(int i = 0; i < 3; i++){
			super.setAtIndex(x, y, z, i, (float) direction.getElement(i));
		}
		
	}
	
	/**
	 * returns the world vector at a given voxel point
	 * @param x
	 * @param y
	 * @param z
	 * @return - vector
	 */
	public SimpleVector getSimpleVectorAtIndex(int x, int y, int z) {

		SimpleVector myVec = new SimpleVector(3);
		
		for (int coord = 0; coord < 3; coord++){ // 3 world coordinates x,y and z
					myVec.setElementValue(coord, getAtIndex(x, y, z, coord));
	}
		return myVec;
	}
	
	
	
	public void writeToVectorField(String pathFiberVTK){

		double max = 0;
		
		ArrayList<SimpleVector> usedPoints = new ArrayList<SimpleVector>();
		
		for (int x = 0; x < imgSizeX; x++){
			for (int y = 0; y < imgSizeY; y++){
				for (int z = 0; z < imgSizeZ; z++){

				
				SimpleVector direction = getSimpleVectorAtIndex(x,y,z);
				double length = direction.normL2();
				
				if(length > max){
					max = length;
				}

				if(length != 0){

				SimpleVector myVec = new SimpleVector(x,y,z,direction.getElement(0),direction.getElement(1),direction.getElement(2),length);
					
				usedPoints.add(myVec);
				
				}
					
				
				}
			}
		}
		
		

		try{
			FileOutputStream foStream = new FileOutputStream(pathFiberVTK);
			BufferedWriter bufWriter = new BufferedWriter(new OutputStreamWriter(foStream));

			bufWriter.write("# vtk DataFile Version 2.0");
			bufWriter.write(System.getProperty( "line.separator" ));
			bufWriter.write("Unstructured Grid Example");
			bufWriter.write(System.getProperty( "line.separator" ));
			bufWriter.write("ASCII");
			bufWriter.write(System.getProperty( "line.separator" ));
			bufWriter.write("DATASET UNSTRUCTURED_GRID");
			bufWriter.write(System.getProperty( "line.separator" ));
			
			// TODO Not sure why that
			bufWriter.write("POINTS " + usedPoints.size() +" float");
			bufWriter.write(System.getProperty( "line.separator" ));

			// First write positions

			for(int i = 0; i < usedPoints.size(); i++){
				SimpleVector curVec = usedPoints.get(i);
				bufWriter.write(curVec.getElement(0) +" " +curVec.getElement(1) +" " +curVec.getElement(2));
				bufWriter.write(System.getProperty( "line.separator" ));
			}
			

			
			bufWriter.write("CELLS " + usedPoints.size() + " " +2*usedPoints.size());
			bufWriter.write(System.getProperty( "line.separator" ));
			
			int indCells = 0;

			for(int i = 0; i < usedPoints.size(); i++){
						bufWriter.write("1 " +indCells);
						bufWriter.write(System.getProperty( "line.separator" ));
						indCells = indCells + 1;
			}
			
			bufWriter.write("CELL_TYPES " + usedPoints.size());
			bufWriter.write(System.getProperty( "line.separator" ));
			
			for(int i = 0; i < usedPoints.size(); i++){
						bufWriter.write("1");
						bufWriter.write(System.getProperty( "line.separator" ));
			}
			
				
				bufWriter.write("POINT_DATA  " + usedPoints.size());
				bufWriter.write(System.getProperty( "line.separator" ));
				bufWriter.write("SCALARS scalarvalue float 1");
				bufWriter.write(System.getProperty( "line.separator" ));
				bufWriter.write("LOOKUP_TABLE default");
				bufWriter.write(System.getProperty( "line.separator" ));
				
				// Write scalars
				for(int i = 0; i < usedPoints.size(); i++){
					bufWriter.write(usedPoints.get(i).getElement(6) +"");
					bufWriter.write(System.getProperty( "line.separator" ));
				}
				
	
				bufWriter.write("VECTORS firstDirection float");
			
				// Write scalars
				for(int i = 0; i < usedPoints.size(); i++){
					SimpleVector normalizedVec = usedPoints.get(i).getSubVec(3, 3);
					normalizedVec.multiplyBy(1/max);
					bufWriter.write(System.getProperty( "line.separator" ));
					bufWriter.write(normalizedVec.getElement(0) + " " + normalizedVec.getElement(1)+" "+"" +normalizedVec.getElement(2));
				}

			bufWriter.flush();
			foStream.close();
			}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
