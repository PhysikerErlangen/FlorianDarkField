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
	SimpleVector[] scatteringVectors;
	
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
		 initializeScatteringVectors();	
		 
		 // Calculate Weight Matrix
		 calculateWeightMatrix();
		 
	}
	
	
	
	
	// Initialize the scattering vectors
	public void initializeScatteringVectors(){
		// For the first test implementations
		if(this.numScatterVectors == 3){
			
			this.scatteringVectors = new SimpleVector[3];
			
			this.scatteringVectors[0] = new SimpleVector(1f,0f,0f); 
			this.scatteringVectors[1] = new SimpleVector(0f,1f,0f);
			this.scatteringVectors[2] = new SimpleVector(0f,0f,1f);
		}	
		
		if(this.numScatterVectors == 7){
			
			this.scatteringVectors = new SimpleVector[7];
			
			this.scatteringVectors[0] = new SimpleVector(1f,0f,0f).normalizedL2(); 
			this.scatteringVectors[1] = new SimpleVector(0f,1f,0f).normalizedL2();
			this.scatteringVectors[2] = new SimpleVector(0f,0f,1f).normalizedL2();
			this.scatteringVectors[3] = new SimpleVector(1f,1f,1f).normalizedL2(); 
			this.scatteringVectors[4] = new SimpleVector(1f,1f,-1f).normalizedL2();
			this.scatteringVectors[5] = new SimpleVector(1f,-1f,1f).normalizedL2();
			this.scatteringVectors[6] = new SimpleVector(-1f,1f,1f).normalizedL2();
			
		}
		
		
		
		else if(this.numScatterVectors == 1){
			
			this.scatteringVectors = new SimpleVector[1];
			
			this.scatteringVectors[0] = new SimpleVector(1,0,1);
			
			 
		}
		
	}

	
	
//	public int getIndex(int uPixel, int vPixel, int proj, int scatterVec){
//		
//		return 0;
//		
//	}
	
	public double getWeight(int projIndex, int scatterVector){
		// TODO DO SOME CALCULATIONS
		return weights[projIndex][scatterVector];
	}

		
	// u,v pixel coordinate, projIndex is projection index, scatterVector current scatter direction
	public void setWeight(int projIndex, int scatterIndex,double val){
		weights[projIndex][scatterIndex] = val;
	}
	

	// Calculate scalar product of two vectors
	public static double scalarProduct(SimpleVector vec1, SimpleVector vec2) {
		assert vec1.getLen() != vec2.getLen() : new IllegalArgumentException(
				"Both vectors have to have same size!");
		double sum = 0;
		for (int k = 0; k < vec1.getLen(); k++) {
			sum += vec1.getElement(k) * vec2.getElement(k);
		}
		return sum;
	}

	// Calculates the cross product between to vectors
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
	// rayDir: Direction of the ray
	// scatterDirection: Direction of the particular Scattering Vector
	// Sensitivity: sensitivity direction of grating interferometer

	public double calculateWeightFactor(SimpleVector rayDir,
			SimpleVector scatterDirection, SimpleVector sensitivity) {
		SimpleVector help = crossProduct(rayDir, scatterDirection);
		double cross = help.normL2();
		double scalar = scalarProduct(scatterDirection, sensitivity);
		double weight = cross * scalar;
		return weight * weight;
	}
	
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

	
	private SimpleVector getScatterVector(int scatterIndex){
		return scatteringVectors[scatterIndex];
	}
	
	private SimpleVector getRayVector(int curTheta){
	
		double theta = deltaTheta * curTheta;

		double cosTheta = Math.cos(theta);
		double sinTheta = Math.sin(theta);
		
		// Second point is just added to p1 but slight perpendicular to it
		SimpleVector rayDir = calculateRotatedVector(-sinTheta,cosTheta,0).getAbstractVector();
		
		return rayDir;
	}
	
	
	private SimpleVector getSensitivityVector(int curTheta){
		
		double theta = deltaTheta * curTheta;

		double cosTheta = Math.cos(theta);		
		double sinTheta = Math.sin(theta);
		
		
		SimpleVector sensitivityVector = calculateRotatedVector(cosTheta,sinTheta,0).getAbstractVector();
		
		return sensitivityVector;
	
}

	public SimpleMatrix getScatterVectorsAsMatrix(){
		
		SimpleMatrix scatterMatrix = new SimpleMatrix(3,scatteringVectors.length );
		
		for(int channel = 0; channel < numScatterVectors; channel++){
			scatterMatrix.setColValue(channel, scatteringVectors[channel]);
		}

		return scatterMatrix;
		
	}
	
	public static void calcScatterDirections(DarkField3DTensorVolume myVolume, SimpleMatrix scatterMatrix, SimpleVector scatterWeights){
		
		for(int x = 0; x <myVolume.imgSizeX; x++){
			for(int y = 0; y <myVolume.imgSizeY; y++){
				for(int z = 0; z <myVolume.imgSizeZ; z++){
					// Initializes the PCA objects
					DarkFieldPCA myPCA = new DarkFieldPCA(scatterMatrix,scatterWeights);
					// Performs PCA
					myPCA.run();
					// Extract Scatter Direction as smallest component of the ellipsoid.
					
				} // End loop z
			} // End loop y
		} // End loop z
		
		
		
	}
	
	
	

}
