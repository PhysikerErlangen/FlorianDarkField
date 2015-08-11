package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;

public class DarkFieldPointCloudViewerTest {

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void test() {
			PointND p1 = new PointND(1,2,3);
			PointND p2 = new PointND(-1,2,-2);
			PointND p3 = new PointND(1,-2,3);
			
			ArrayList<PointND> points = new ArrayList<PointND>();
			points.add(p1);
			points.add(p2);
			points.add(p3);
			
			DarkFieldPointCloudViewer myMaker = new DarkFieldPointCloudViewer(points);

	}

}
