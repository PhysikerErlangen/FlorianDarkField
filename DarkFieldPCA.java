// This code was developed in a collaboration with ECAP, Erlangen, Germany.
// This part of the code is not to be published under GPL before Oct 31st 2017.
// author@ Florian Schiffers August 3st, 2015
//

package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import ij.ImageJ;
import ij.gui.Plot;

import java.util.ArrayList;

import edu.mines.jtk.dsp.Eigen;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.geometry.shapes.activeshapemodels.PCA;
import edu.stanford.rsl.conrad.geometry.shapes.mesh.DataMatrix;
import edu.stanford.rsl.conrad.numerics.DecompositionSVD;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;

public class DarkFieldPCA{

	
	// Test implemtation of the PCA Method
	public static void main(String[] args){
		
		double[][] multi = new double[][]{
				  { 1, 1, 0},
				  { 1, 1, 0},
				  { 0, 2, 3},
				  { 1, 1, 0},
				  { 1, 1, 1},
				  { 5, 4, 4},
				  { 5, 4, 4},
				}; 
	
		SimpleVector myScatterWeights = new SimpleVector(1,1,1,1,1,1,1);
		
		SimpleMatrix scatterDirections = new SimpleMatrix(multi);
		scatterDirections = scatterDirections.transposed();
	
		DarkFieldPCA myPCA = new DarkFieldPCA(scatterDirections, myScatterWeights);
		
		myPCA.run();
		
		myPCA.showDataPoints();
		
	}
	
	// Contains the different scatterDirections (is constant and given as input)
	SimpleMatrix scatterDirections;
	
	// Scatter points contains the scatter directions of the respective 
	// Voxel which are already weighted by the earlier calculated
	// scatter coefficients. These points are then used for calculating
	// the principal axes by PCA
	SimpleMatrix scatterPoints;
	
	
	/**
	 * Data-Set consists of a set of many vectors of the same dimensions. The number of points is saved in numPoints
	 */
	private int numPoints;
	
	
	/**
	 * Data-Set consists of a set of many vectors of the same dimensions. This dimension is saved in dimension
	 */
	private int dimension;
	
	/**
	 * Matrix containing the Eigenvectors of the covariance matrix after singular value decomposition.
	 */
	private SimpleMatrix eigenVectors;
			
	
	/**
	 * @return the eigenVectors
	 */
	public SimpleMatrix getEigenVectors() {
		return eigenVectors;
	}

	/**
	 * contains the calculated scatter coefficient by tensor reconstruction
	 */
	private SimpleVector myScatterWeights;
	
	/**
	 * covMatrix contains the entries covariance matrix of the given the data set
	 */
	private SimpleMatrix covMatrix;
	
	
	
	/**
	 * Array containing the eigenvalues of the covariance matrix after singular value decomposition.
	 */
	private SimpleVector eigenValues;
	
	
	/**
	 * @return the eigenValues
	 */
	public SimpleVector getEigenValues() {
		return eigenValues;
	}

	/**
	 * Constructor for DarkFieldPCA.
	 * CAUTION: DarkFieldPCA needs special data, as no mean shift is performed, as the
	 * data points are already assumed to have a mean of 0.
	 */
	public DarkFieldPCA (SimpleMatrix scatterDirections, SimpleVector myScatterWeights){

		assert(scatterDirections.getCols() != myScatterWeights.getLen()) : new Exception("Dimensions do not match in DarkFieldPCA.");
		
		// Save reference to the predefined scatterDirections
		this.scatterDirections = scatterDirections;
		
		this.myScatterWeights = myScatterWeights;
		
		// Get Dimensionality of one scatter direction
		dimension = scatterDirections.getRows();
		
		// Get number of scatter points
		// 2*numSamples because we take + and - direction of the scatter direction
		// for calculating the PCA due to stability reasons (see paper of Vogel)
		numPoints = 2*myScatterWeights.getLen();
		
		// Calc the matrix used for PCA with dimension: 3 X 2*Num_Points
		
		
		
		
	}
		
	// Calculates the different scatter directions with their corresponding weights. 
	private void calcSetOfScatterPoints(){
		
		scatterPoints = new SimpleMatrix(dimension,numPoints);
		
		for (int direction = 0; direction < myScatterWeights.getLen(); direction++){
			
			double myWeight = myScatterWeights.getElement(direction);
			SimpleVector vec1 = scatterDirections.getCol(direction);
			vec1.multiplyBy(myWeight);
			SimpleVector vec2 = vec1.multipliedBy(-1f);
			
			scatterPoints.setColValue(2*direction, vec1);
			scatterPoints.setColValue(2*direction + 1, vec2);
		}
	} // END calcSetOFScatterPoints
	
	
	/**
	 * Calculates covariance matrix
	 */
	
	private void calcCovarianceMatrix(){
		
		assert(scatterPoints != null) : new Exception("Initialize data array fist.");
		
		// Init covariance matrix
		covMatrix = new SimpleMatrix(dimension,dimension);

		// Loop through all data points to calculate covariance matrix
		for(int ind = 0; ind < numPoints; ind++){
		// Get current vector
		SimpleVector curVec = scatterPoints.getCol(ind);
		// Calculate out product
		SimpleMatrix outerMatrix = SimpleOperators.multiplyOuterProd(curVec, curVec);
		// Add current matrix to covariance matrix
		covMatrix.add(outerMatrix);
		}
		
		covMatrix.divideBy((numPoints-1));
		
	}
	
	/**
	 * Plots the data in the array over its array index.
	 * Written my Matthias Unberath. Adapted by Florian Schiffers
	 * @param data The data to be plotted.
	 */
	private void plot(double[] data){
		
			new ImageJ();
			Plot plot = VisualizationUtil.createPlot(data, "Singular values of data matrix", "Singular value", "Magnitude");
			plot.show();
		
	}
	
	
	public void showDataPoints(){
		
		DarkFieldPointCloudViewer myViewer = new DarkFieldPointCloudViewer(scatterPoints);
		
		myViewer.showPoints();
		
	}
	
	
	
	/**
	 * Performs the principal component analysis on the data-set.
	 */
	public void run(){
		assert(scatterPoints != null) : new Exception("Initialize data array fist.");
		
		// First calculate the "Dataset" that is used for PCA. This dataset has already by definition mean 0
		calcSetOfScatterPoints();
		
		// Calculate covariance matrix of given dataset
		calcCovarianceMatrix();
		
		DecompositionSVD svd = new DecompositionSVD(covMatrix);
		
		//
		// plot(svd.getSingularValues());
		
		eigenValues  =  new SimpleVector( svd.getSingularValues());
		eigenVectors =  svd.getU();
	
	}
	
	
}
