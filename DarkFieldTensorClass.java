package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;

import edu.stanford.rsl.conrad.geometry.General;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class DarkFieldTensorClass{

	// Defines dimension of volume box
	protected int imgSizeX;
	protected int imgSizeY;
	protected int imgSizeZ;
	
	private ArrayList<DarkFieldVectorField> eigenVectorList;

	/**
	 * @param imgSizeX
	 * @param imgSizeY
	 * @param imgSizeZ
	 * @param spacing_world
	 */
	public DarkFieldTensorClass(int imgSizeX, int imgSizeY,int imgSizeZ, double[] spacing_world, double[] origin_world){
		
		this.imgSizeX = imgSizeX;
		this.imgSizeY = imgSizeY;
		this.imgSizeZ = imgSizeZ;

		// contain the eigenvectors (4th column is always the respective eigenvalue)
		DarkFieldVectorField eigenVec1 = new DarkFieldVectorField(imgSizeX, imgSizeY, imgSizeZ, spacing_world, origin_world);
		DarkFieldVectorField eigenVec2 = new DarkFieldVectorField(imgSizeX, imgSizeY, imgSizeZ, spacing_world, origin_world);
		DarkFieldVectorField eigenVec3 = new DarkFieldVectorField(imgSizeX, imgSizeY, imgSizeZ, spacing_world, origin_world);
		
		eigenVectorList = new ArrayList<DarkFieldVectorField>(3);
		eigenVectorList.add(eigenVec1);
		eigenVectorList.add(eigenVec2);
		eigenVectorList.add(eigenVec3);
		

	}
	
	public DarkFieldVectorField getFiberDirection(){
		return eigenVectorList.get(2);
	}
	
	public void setData(int x, int y, int z, DarkFieldPCA myPCA){
		
		
		SimpleMatrix eigenVectors = myPCA.getEigenVectors();
		SimpleVector eigenValues = myPCA.getEigenValues();
		
		for(int i = 0; i < eigenVectorList.size(); i++){
			SimpleVector eigenVec = eigenVectors.getCol(i);
			double eigenVal = eigenValues.getElement(i);
			eigenVec = eigenVec.multipliedBy(eigenVal);
			eigenVectorList.get(i).setVector(x, y, z, eigenVec);
		}
		

//		/**
//		 * Threshold that checks, if 3 component of eigenvalues is too small
//		 * If 3 component is too small, don't consider it as a fiber orientation
//		 * and ignore it
//		 */
//		double th = 1E-10;

//		SimpleVector fiberDir;
//		if(myPCA.getEigenValues().getElement(2)<th){
//			fiberDir = new SimpleVector(3);
//		}else{
//			fiberDir = myPCA.getEigenVectors().getCol(2).normalizedL2();
//			fiberDir.multiplyBy(myPCA.getEigenValues().getElement(2));
//		}
		
		
		
	}
	
	
	/**
	 * Helper class, containing 3 index points representing
	 * the x,y,z coordinates of a given mesh point
	 */
	class Index3D {

		int x;
	    int y;
	    int z;
		
		Index3D(int x, int y, int z){
			
			this.x = x;
			this.y = y;
			this.z = z;
		}
		}
	
	
	/**
	 * @param pathFiberVTK
	 */
	public void writeToVectorField(String pathFiberVTK){

		
		ArrayList<Index3D> indexListe = new ArrayList<Index3D>();
		
		for (int x = 0; x < imgSizeX; x++){
			for (int y = 0; y < imgSizeY; y++){
				for (int z = 0; z < imgSizeZ; z++){
				
				// Get EigenVec of smallest eigenValue
				SimpleVector direction = eigenVectorList.get(2).getSimpleVectorAtIndex(x,y,z);
				
				double length = direction.normL2();
				
				if(length != 0){
					Index3D curIndex = new Index3D(x, y, z);
					indexListe.add(curIndex);
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
			bufWriter.write("POINTS " + indexListe.size() +" float");
			bufWriter.write(System.getProperty( "line.separator" ));

			// First write positions

			for(int poinIdx = 0; poinIdx < indexListe.size(); poinIdx++){
				Index3D index3D = indexListe.get(poinIdx);
				bufWriter.write(index3D.x +" " +index3D.y +" " +index3D.z);
				bufWriter.write(System.getProperty( "line.separator" ));
			}
			

			
			bufWriter.write("CELLS " + indexListe.size() + " " +2*indexListe.size());
			bufWriter.write(System.getProperty( "line.separator" ));
			
			int indCells = 0;

			for(int poinIdx = 0; poinIdx < indexListe.size(); poinIdx++){
						bufWriter.write("1 " +indCells);
						bufWriter.write(System.getProperty( "line.separator" ));
						indCells = indCells + 1;
			}
			
			bufWriter.write("CELL_TYPES " + indexListe.size());
			bufWriter.write(System.getProperty( "line.separator" ));
			
			for(int pointIdx = 0; pointIdx < indexListe.size(); pointIdx++){
						bufWriter.write("1");
						bufWriter.write(System.getProperty( "line.separator" ));
			}
			
				
				bufWriter.write("POINT_DATA  " + indexListe.size());
				bufWriter.write(System.getProperty( "line.separator" ));
				bufWriter.write("SCALARS eigenValues float 3");
				bufWriter.write(System.getProperty( "line.separator" ));
				bufWriter.write("LOOKUP_TABLE default");
				bufWriter.write(System.getProperty( "line.separator" ));
				
				// Write scalars
				for(int pointIdx = 0; pointIdx < indexListe.size(); pointIdx++){
					
					Index3D coordIdx = indexListe.get(pointIdx);
					SimpleVector vec1 = eigenVectorList.get(0).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
					SimpleVector vec2 = eigenVectorList.get(1).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
					SimpleVector vec3 = eigenVectorList.get(2).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
					
					bufWriter.write(vec1.normL2() +" " + vec2.normL2() +" " + vec3.normL2());
					bufWriter.write(System.getProperty( "line.separator" ));
				}
				
				bufWriter.write("SCALARS scalarDirection float 1");
				bufWriter.write(System.getProperty( "line.separator" ));
				bufWriter.write("LOOKUP_TABLE default");
				bufWriter.write(System.getProperty( "line.separator" ));
				
				// Write scalars
				for(int pointIdx = 0; pointIdx < indexListe.size(); pointIdx++){

					Index3D coordIdx = indexListe.get(pointIdx);
					SimpleVector fiberDir = eigenVectorList.get(2).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
					// Set z component to 0, as we don't care about this in the visualization.
					fiberDir.setElementValue(2, 0);
					
					/*
					 * Reference vector for scalar value for plot visualization in vtk.
					 * Scalar value is calculated by scalar product with two vectors 
					 */
					SimpleVector refVec = new SimpleVector(1f,0,0);
					double angle = Math.abs(General.angle(fiberDir,refVec));
					
					angle = Math.abs(angle);
					if(angle >= Math.PI){
						
					}
					
					bufWriter.write(angle+"");
					bufWriter.write(System.getProperty( "line.separator" ));
					
				}
				
	
				
			
				for(int eigIdx = 0; eigIdx < 3; eigIdx++){
				
				bufWriter.write("VECTORS eigenVector" +eigIdx  +" float");
					
				// Write scalars
				for(int pointIdx = 0; pointIdx < indexListe.size(); pointIdx++){
			
					Index3D coordIdx = indexListe.get(pointIdx);
					SimpleVector curEigVec = eigenVectorList.get(eigIdx).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
					bufWriter.write(System.getProperty( "line.separator" ));
					bufWriter.write(curEigVec.getElement(0) + " " + curEigVec.getElement(1)+" "+"" +curEigVec.getElement(2));
				}
				bufWriter.write(System.getProperty( "line.separator" ));
				}
				
				
					
				
					bufWriter.write("TENSORS tensorData float");
						
					// Write scalars
					for(int pointIdx = 0; pointIdx < indexListe.size(); pointIdx++){
				
						Index3D coordIdx = indexListe.get(pointIdx);
						
						SimpleVector eigVec1 = eigenVectorList.get(0).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
						SimpleVector eigVec2 = eigenVectorList.get(1).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
						SimpleVector eigVec3 = eigenVectorList.get(2).getSimpleVectorAtIndex(coordIdx.x, coordIdx.y, coordIdx.z);
						
						
						SimpleMatrix mat1 = SimpleOperators.multiplyOuterProd(eigVec1, eigVec1).multipliedBy(1/eigVec1.normL2());
						SimpleMatrix mat2 = SimpleOperators.multiplyOuterProd(eigVec2, eigVec2).multipliedBy(1/eigVec2.normL2());
						SimpleMatrix mat3 = SimpleOperators.multiplyOuterProd(eigVec3, eigVec3).multipliedBy(1/eigVec3.normL2());
						
						SimpleMatrix cov = SimpleOperators.add(mat1,mat2,mat3);
						
						bufWriter.write(System.getProperty( "line.separator" ));
						bufWriter.write(cov.getElement(0, 0) + " " + cov.getElement(1, 0) + " " + cov.getElement(2, 0) + " ");
						bufWriter.write(System.getProperty( "line.separator" ));
						bufWriter.write(cov.getElement(0, 1) + " " + cov.getElement(1, 1) + " " + cov.getElement(2, 1) + " ");
						bufWriter.write(System.getProperty( "line.separator" ));
						bufWriter.write(cov.getElement(0, 2) + " " + cov.getElement(1, 2) + " " + cov.getElement(2, 2) + " ");
						bufWriter.write(System.getProperty( "line.separator" ));
						
//						bufWriter.write(System.getProperty( "line.separator" ));
//						bufWriter.write(eigVec1.getElement(0) + " " + eigVec1.getElement(1) + " " + eigVec1.getElement(2));
//						bufWriter.write(System.getProperty( "line.separator" ));
//						bufWriter.write(eigVec2.getElement(0) + " " + eigVec2.getElement(1) + " " + eigVec2.getElement(2));
//						bufWriter.write(System.getProperty( "line.separator" ));
//						bufWriter.write(eigVec3.getElement(0) + " " + eigVec3.getElement(1) + " " + eigVec3.getElement(2));
//						bufWriter.write(System.getProperty( "line.separator" ));
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


 
