package edu.stanford.rsl.science.darkfield.FlorianDarkField;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

public class DarkField3DTensorVolumeTest {

	DarkField3DTensorVolume volA,volB;
	
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		
		int imgSizeX = 1;
		int imgSizeY = 1;
		int imgSizeZ = 1;
		
		int numChannels = 2;
		
		double[] spacing_world = {1,1,1};
		double[] origin_world = {-1,-1,-1};
		
		volA = new DarkField3DTensorVolume(imgSizeX, imgSizeY, imgSizeZ, numChannels, spacing_world, origin_world);
		volB = new DarkField3DTensorVolume(imgSizeX, imgSizeY, imgSizeZ, numChannels, spacing_world, origin_world);
		
		volA.setAtIndex(0, 0, 0, 0, 1);
		volA.setAtIndex(0, 0, 0, 1, 2);
		
		volB.setAtIndex(0, 0, 0, 0, 2);
		volB.setAtIndex(0, 0, 0, 1, 4);
		
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testSub() {
		
		DarkField3DTensorVolume sub = DarkField3DTensorVolume.sub(volB, volA);
		
		assertTrue("Sub does not work fine.",sub.getAtIndex(0, 0, 0, 0) == 1 );
		assertTrue("Sub does not work fine.",sub.getAtIndex(0, 0, 0, 1) == 2);
		
	}
	
	@Test
	public void testSubReference(){
		
		volA.sub(volA);
		
		assertTrue("Sub does not work fine.",volA.getAtIndex(0, 0, 0, 0) == 0 );
		assertTrue("Sub does not work fine.",volA.getAtIndex(0, 0, 0, 1) == 0 );
		
		
	}
	

}
