package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import java.util.ArrayList;

import edu.stanford.rsl.apps.gui.opengl.PointCloudViewer;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;

public class DarkFieldPointCloudViewer {

	ArrayList<PointND> points;
	
	PointCloudViewer pointCloudViewer;
	
	public static void main(String[] args){
		
	
		
		PointND p1 = new PointND(1,2,3);
		PointND p2 = new PointND(-1,2,-2);
		PointND p3 = new PointND(1,-2,3);
		
		ArrayList<PointND> points = new ArrayList<PointND>();
		points.add(p1);
		points.add(p2);
		points.add(p3);
		
		DarkFieldPointCloudViewer myMaker = new DarkFieldPointCloudViewer(points);
		
	}
	
	
	public void showPoints(){
		showPoints(true);
	}
	
	public void showPoints(boolean flag){
		pointCloudViewer.setVisible(flag);
	}
	
	
	public DarkFieldPointCloudViewer(SimpleMatrix myPoints) {
	this(myPoints,"");
	}
	
	public DarkFieldPointCloudViewer(SimpleMatrix myPoints, String title) {
	   
		assert(myPoints.getRows()==3) : new Exception("Dimension of data points has to be 3!");
		assert(myPoints.getCols()!=0) : new Exception("Dimension of data points has to be larger 0!");
	   
		ArrayList<PointND> myPointList = new ArrayList<PointND>();
		
		for(int pointInd = 0; pointInd < myPoints.getCols(); pointInd++){
			PointND myPoint = new PointND(myPoints.getCol(pointInd));
			myPointList.add(myPoint);
		}
		
		initData(myPointList,title);
	}

	
public void initData(ArrayList<PointND> points, String title){
	
	this.points = points;
	// Open new PointCloudViewer Object
	pointCloudViewer = new PointCloudViewer(title, points);
}
	
	public DarkFieldPointCloudViewer(ArrayList<PointND> points){
	this(points,"");
	}
	
	public DarkFieldPointCloudViewer(ArrayList<PointND> points, String title){
	initData(points, title);	
	}
	
	
	
}
