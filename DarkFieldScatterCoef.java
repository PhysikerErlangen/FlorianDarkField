package edu.stanford.rsl.science.darkfield.FlorianDarkField;


import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.science.overexposure.CrossCalibration;



public class DarkFieldScatterCoef extends  DarkFieldTensorGeometry{

	public static void main(String[] args){
		
		SimpleVector vec1 = new SimpleVector(1,0,0);
		SimpleVector vec1_2 = new SimpleVector(2,0,0);
		SimpleVector vec2 = new SimpleVector(0,1,0);
		SimpleVector vec3 = new SimpleVector(0,0,1);
		
		SimpleVector vec4 = new SimpleVector(0,1,1);
		SimpleVector vec5 = new SimpleVector(1,1,1);
		
		SimpleVector vec6 = new SimpleVector(-2.3,1.3,7.3); 
		SimpleVector vec7 = new SimpleVector(1.1,8.3,-11);
		
		SimpleVector res;
		
		// Check if cross product is working correctly
		res = crossProduct(vec1, vec2);
		System.out.println(vec1 + "x" +vec2 +" =  " +res);
		
		res = crossProduct(vec1, vec3);
		System.out.println(res);
		
		res = crossProduct(vec1_2, vec1);
		System.out.println(res);
		
		res = crossProduct(vec4, vec5);
		System.out.println(res);
		
		res = crossProduct(vec6, vec7);
		System.out.println(res);
		 
		double resD;
		// Check if scalar product is working nice
		// Check if cross product is working correctly
		resD = scalarProduct(vec1, vec2);
		System.out.println(vec1 + "*" +vec2 +" =  " +resD);
		
		resD = scalarProduct(vec1, vec3);
		System.out.println(resD);
		
		resD = scalarProduct(vec1_2, vec1);
		System.out.println(resD);
		
		resD = scalarProduct(vec4, vec5);
		System.out.println(resD);
		
		resD = scalarProduct(vec6, vec7);
		System.out.println(resD);
		
		
	}
	
	// Declares array of weighting factors
	private double[][] weights;

	int numElements;

	int numScatterVectors;
	
	// Number of Scattering vectors
	SimpleMatrix scatteringVectors;
	
	/**
	 * @param conf
	 * @param numScatterVectors
	 */
	public DarkFieldScatterCoef( Configuration conf, int numScatterVectors) {

		// Call super constructor which contains whole geometry etc.
		super(conf,numScatterVectors);
		
		// Number of weights: NumberPixels * NumberProjections *
		// NumberScatterVectors
		// int numWeights = maxV * maxU * maxProjs * numScatterVectors;
		// Allocate data for the weights data (might be really, really huge)
		
		System.out.println("Start allocating scatter weights.");
		
		this.numScatterVectors = numScatterVectors;
		
		// numElements = maxProjs * numScatterVectors;
		
		 // weights = new double[maxU_index][maxV_index][maxProjs][numScatterVectors];
		 weights = new double[maxTheta_index][numScatterVectors];
		 
		 System.out.println("Scatter weights allocated. Total number of elements: " +(maxTheta_index*numScatterVectors));


		 // initialize the scattering vectors of the grid
		 scatteringVectors = DarkFieldScatterDirection.getScatterDirectionMatrix(numScatterVectors);	
		 
		 // Calculate Weight Matrix
		 calculateWeightMatrix();
		 
	}
	
	
	
	/**
	 * @param curTheta
	 * @param scatterVector
	 * @return
	 */
	public double getWeight(int curTheta, int scatterVector){
		// TODO DO SOME CALCULATIONS
		return weights[curTheta][scatterVector];
	}


	/**
	 * @param curTheta
	 * @param scatterIndex
	 * @param val
	 */
	public void setWeight(int curTheta, int scatterIndex,double val){
		weights[curTheta][scatterIndex] = val;
	}
	
 
	/**
	 * Calculate scalar product of two vectors
	 * @param vec1
	 * @param vec2
	 * @return
	 */
	public static double scalarProduct(SimpleVector vec1, SimpleVector vec2) {
		assert vec1.getLen() != vec2.getLen() : new IllegalArgumentException(
				"Both vectors have to have same size!");
		double sum = 0;
		for (int k = 0; k < vec1.getLen(); k++) {
			sum += vec1.getElement(k) * vec2.getElement(k);
		}
		return sum;
	}

 
	/**
	 * Calculates the cross product between to vectors
	 * @param vec1
	 * @param vec2
	 * @return
	 */
	public static SimpleVector crossProduct(SimpleVector vec1, SimpleVector vec2) {
		// Assert if vector length is not the same
		assert (vec1.getLen() != 3) || (vec2.getLen() != 3) : new IllegalArgumentException(
				"Length of both vectors has to be 3");

		double newA = vec1.getElement(1) * vec2.getElement(2)
				- vec1.getElement(2) * vec2.getElement(1);
		double newB = vec1.getElement(2) * vec2.getElement(0)
				- vec1.getElement(0) * vec2.getElement(2);
		double newC = vec1.getElement(0) * vec2.getElement(1)
				- vec1.getElement(1) * vec2.getElement(0);

		return new SimpleVector(newA, newB, newC);

	}

	// calculateWeightFactor calculates the weight factor which are defined in
	// maleckis paper
	// rayDir: 
	// scatterDirection: 
	// Sensitivity: 

	/**
	 * 
	 * @param rayDir - Direction of the ray
	 * @param scatterDirection - Direction of the particular Scattering Vector
	 * @param sensitivity - sensitivity direction of grating interferometer
	 * @return
	 */
	public double calculateWeightFactor(SimpleVector rayDir,
			SimpleVector scatterDirection, SimpleVector sensitivity) {
		SimpleVector help = crossProduct(rayDir, scatterDirection);
		double cross = help.normL2();
		double scalar = scalarProduct(scatterDirection, sensitivity);
		double weight = cross * scalar;
		return weight * weight;
	}
	
	/**
	 * Calculates the whole weight matrix according to Eq. (2)
	 * of Paper: Vogel - Constrained X-Ray tensor tomography reconstruction
	 * Here parallel Beam structure is assumed, which is why the actual
	 * dimension of 
	 */
	public void calculateWeightMatrix() {

		/* Assume that object is static and Detector is moved around object
		 * 
		 * scatterVectors are constant and don't change direction
		 * Ray direction has to be calculated new for every projection
		 * sensitivity is constant
		 * 
		 */
				
		// Loop through all projections
		for (int curTheta = 0; curTheta < this.maxTheta_index; curTheta++) {
			// Sensitivity is dependent on the projection
			SimpleVector sensitivity = this.getSensitivityVector(curTheta);
		
		// Loop through all scatter vectors
		for (int scatterIndex = 0; scatterIndex < numScatterVectors; scatterIndex++) {
			SimpleVector scatterDirection = getScatterVector(scatterIndex);

				
						// Direction of the ray is dependent on projection
						SimpleVector rayDir = getRayVector(curTheta); 
						// calculate current weight factor
						double curWeight = calculateWeightFactor(rayDir,scatterDirection, sensitivity);
						setWeight(curTheta,scatterIndex,curWeight);
				
						// System.out.println("Scatter Vector " +scatterIndex +" theta: " +curTheta + " curWeight: " +curWeight);
						
			}
		}
	}

	
	/**
	 * returns scatter sample direction of the particular channel
	 * @param scatterIndex
	 * @return
	 */
	private SimpleVector getScatterVector(int scatterIndex){
		return scatteringVectors.getCol(scatterIndex);
	}
	
	/**
	 * Return the x-Ray direction at a projection angle theta
	 * @param curTheta
	 * @return
	 */
	private SimpleVector getRayVector(int curTheta){
	
		double theta = deltaTheta * curTheta;

		double cosTheta = Math.cos(theta);
		double sinTheta = Math.sin(theta);
		
		// Second point is just added to p1 but slight perpendicular to it
		SimpleVector rayDir = calculateRotatedVector(-sinTheta,cosTheta,0).getAbstractVector();
		
		return rayDir;
	}
	
	
	/**
	 * Returns the sensitivity vector at a projection angle theta
	 * @param curTheta
	 * @return
	 */
	private SimpleVector getSensitivityVector(int curTheta){
		
		double theta = deltaTheta * curTheta;

		double cosTheta = Math.cos(theta);		
		double sinTheta = Math.sin(theta);
		
		
		SimpleVector sensitivityVector = calculateRotatedVector(cosTheta,sinTheta,0).getAbstractVector();
		
		return sensitivityVector;
	
}


	
	
	

}
