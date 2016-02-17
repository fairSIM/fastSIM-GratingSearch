/*
This file is part of fastSIM Grating Search, ported to Java.

This code is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

The code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the code.  If not, see <http://www.gnu.org/licenses/>
*/

package de.bio_photonics.gratingsearch;

import java.util.ArrayList;
import java.util.List;

import org.fairsim.fiji.ImageVector;
import org.fairsim.linalg.Cplx;
import org.fairsim.linalg.Vec2d;
import org.fairsim.linalg.Transforms;
import org.fairsim.sim_algorithm.SimUtils;

import org.fairsim.fiji.DisplayWrapper;

import org.fairsim.utils.ImageDisplay;
import org.fairsim.utils.Tool; 

import ij.ImagePlus;
import ij.ImageStack;

/** Java port of fastSIM SLM grating search algorithm.
The original MATLAB code can be found here, please cite their
publication if you use this software to create SIM gratings:

https://github.com/nanoimaging/fastSIM_GratingSearchforSLM
*/
public class Grating_Search implements ij.plugin.PlugIn {

    // search value space
    final int axmax = 28, aymax = 28;
    final int bxmin = 2,  bymin = 2;
    final int bxmax = 28, bymax = 28;

    // fourier space pxl
    final int fsPxl = 512;

    // SLM #pxl
    int slmPxlX = 1280, slmPxlY = 1024;
    
    
    /** Calculated possible gratings by iterating parameter space.
     *  Each parameter set and resulting grating is tested to meet
     *  conditions: (a) grating period within given limits, (b)
     *  grating is shiftable by a given number of equi-distant phase steps.
     *  
     *  @param gratPerMin Minimal grating period, in #pixels
     *  @param gratPerMax Maximal grating period, in #pixels
     *  @param phaseSteps Number of phase steps
     *  @return List of gratings satisfying the conditions
     *  */
    public List<Grating> calcGrat(
	final double gratPerMin, 
	final double gratPerMax,
	final int phaseSteps ) {

	// List of grating candidates
	List<Grating> candidates = new ArrayList<Grating>(1000);
	
	// outer loop set: loop grating sizes
	for (int ax = 0; ax < axmax; ax ++ ) {
	    System.out.println(String.format(
		"DONE: %d / %d, #candidates %d",ax,axmax, 
		    candidates.size()));
	    for( int ay = -aymax; ay < aymax; ay++ ) {

		// inner loop set: loop grating shifts
		for ( int bx =  bxmin; bx < bxmax; bx++ )
		for ( int by = -bxmax; by < bymax; by++ ) {
		     
		    if (Math.abs(by)<bymin) continue;

		    // create the grating, auto-computes parameters
		    Grating current = new Grating(ax,ay,bx,by);
	   
		    // check if the grating period is within bounds
		    if ( current.gratPer < gratPerMin || 
			 current.gratPer > gratPerMax )
			continue;

		    // check if n equid. pha. can be shifted horizontal
		    if ( current.testStepHorizontal(phaseSteps) ) {
			current.shiftDir = 0;
			fuzzyAdd( current, candidates );
		    }
		    // check if n equid. pha. can be shifted vertical
		    if ( current.testStepVertical(phaseSteps) ) {
			current.shiftDir = 1;
			fuzzyAdd( current, candidates );
		    }


		}
	    }
	}

	return candidates;
    }

    /** Add a grating only if no grating with similar direction
     *  is already present */
    private static boolean fuzzyAdd( Grating a , List<Grating> candidates ) {
	
	// Tolerance when looking for equal direction, in rad
	final double tolDir = 1e-5;
	
	for ( Grating b : candidates )
	    if ( Math.abs( b.gratDir  - a.gratDir ) < tolDir )
		return false;
	
	candidates.add( a );
	return true;
    }


    /** Select combinations of grating orientations.
     *  From a list of possible gratings, generate subsets
     *  with equidistant angles, i.e. for n directions, angle
     *  between gratings is approx. pi/n.
     *
     *  @param candidates   List of grating candidates
     *  @param nrDirs	    Number of directions (typ. 3, 5, ...)
     *  @param tolDirection Allowed deviation from ideal angle, in rad
     *	@return List of matching grating combinations
     *	*/
    public List<Grating []> selectDirs( 
	List<Grating> candidates, 
	final int nrDirs, final double tolDirection ) {

	double sector = Math.PI/nrDirs;

	List<Grating []> ret = new ArrayList<Grating []>();

	// outer loop, select only dirs in the first sector
	for ( Grating c : candidates ) {
	    if (!c.testDirection(0, sector)) continue;

	    // create the array
	    Grating [] col = new Grating[nrDirs];
	    col[0] = c;

	    // see if the other directions can be found within limits
	    for (int i=1; i<nrDirs; i++) {
		double minDist=tolDirection;
		for ( Grating d : candidates ) {
		    double dist =  Math.abs(
			    Grating.wrapPi(d.gratDir) - sector*i - 
			    Grating.wrapPi(c.gratDir) );
		    if ( dist < minDist ) {
			col[i]  = d;
			minDist = dist;
		    }
		}
	    }

	    // check if a complete set was found
	    boolean complete = true;
	    for ( Grating e : col )
		if ( e == null ) complete = false;

	    if (complete)
		ret.add( col );

	}

	return ret;
    }

    /** For a combination of gratings, check if a mask could block unwanted orders. 
     *	@param candidates   List of gratings to check
     *	@param illum	    Illumination vector to apply
     *	@param maxUnwanted  Maximum ratio of amplitude passing through mask
     *	@param maskSize	    Size of mask, in pixel
     *	@param	spatial	    Receives spatial images of illuminated SLM, may be null
     *	@param  fourier	    Receives Fourier images of illuminated SLM, may be null
     * */
    public boolean fourierCheck( Grating [] candidates, Vec2d.Real illum, 
	double maxUnwated, int maskSize,
	ImageDisplay spatial, ImageDisplay fourier ) {

	Vec2d.Cplx     sumFreq = Vec2d.createCplx( fsPxl, fsPxl );
	Vec2d.Cplx [] gratFreq = Vec2d.createArrayCplx( candidates.length, fsPxl, fsPxl );

	Tool.Timer t1 = Tool.getTimer();
	Tool.Timer t2 = Tool.getTimer();
	Tool.Timer t3 = Tool.getTimer();
	Tool.Timer t4 = Tool.getTimer();
	Tool.Timer t5 = Tool.getTimer();

	// Compute the gratings Fourier space and sum them up
	for ( int i=0; i<candidates.length; i++ ) {
	    t1.start();
	    candidates[i].writeToVector( gratFreq[i] );	// get grating as vector 
	    t1.stop();

	    t2.start();
	    gratFreq[i].times( illum );
	    t2.stop();

	    t3.start();
	    gratFreq[i].fft2d(false);			// transform
	    t3.stop();

	    t4.start();
	    Transforms.swapQuadrant(gratFreq[i]);	// quadrant-swap FFT result 
	    sumFreq.add( gratFreq[i] );			// add to full vector
	    t4.stop();

	}

	// store spectrum (if display is set != null)
	t5.start();
	if ( fourier!=null ) {
	    
	    // markers for positions
	    ImageDisplay.Marker [] maskRings = new ImageDisplay.Marker[candidates.length*4];
	    for ( int i=0; i<candidates.length; i++ ) {
		double [] mPos = candidates[i].peakPos(fsPxl);
		maskRings[i*4+0] = new ImageDisplay.Marker(256+mPos[0], 256+mPos[1], 10, 10, true);
		maskRings[i*4+1] = new ImageDisplay.Marker(256-mPos[0], 256-mPos[1], 10, 10, true);
		maskRings[i*4+2] = new ImageDisplay.Marker(256+mPos[0], 256+mPos[1], maskSize, maskSize, false);
		maskRings[i*4+3] = new ImageDisplay.Marker(256-mPos[0], 256-mPos[1], maskSize, maskSize, false);
	    }
	    fourier.addImage(magnitude( sumFreq ), "Spectrum", maskRings);
	}
	t5.stop();

	Tool.trace("Timing: "+t1+t2+t3+t4+t5);
	return true;

    }


    /** Create a Gaussian illumination profile. */
    public Vec2d.Real createGaussIllum( double fwhm, int size) {

	final double sigma = fwhm/2.355;
	Vec2d.Real vec = Vec2d.createReal(size,size);

	for (int y=0; y<size; y++)
	for (int x=0; x<size; x++) {

	    double dist = Math.hypot( y-size/2. , x-size/2. );
	    double fac = Math.exp( -(dist*dist) / (2 * sigma*sigma));
	    vec.set(x,y, (float)fac);
	}
	return vec;
    }

    /** Convolve a Fourier spectrum with the residual SLM structure */
    public void convSLMstructure(Vec2d.Cplx spec) {

	spec.fft2d(true);
	Vec2d.Cplx residual = Vec2d.createCplx(spec);

	final int w = spec.vectorWidth();
	final int h = spec.vectorWidth();

	// structure directly from the original matlab script
	residual.set( 0,   h/2-1, Cplx.Float.one().mult(0.08f));   // top
	residual.set( w-1, h/2-1, Cplx.Float.one().mult(0.08f));   // bottom
	
	residual.set( w/4-1,   0, Cplx.Float.one().mult(0.03f));   // top left
	residual.set( w*3/4-1, 0, Cplx.Float.one().mult(0.03f));   // bottom left
	
	residual.set( w/4-1,   h/2-1, Cplx.Float.one().mult(0.03f));   // top left
	residual.set( w*3/4-1, h/2-1, Cplx.Float.one().mult(0.03f));   // bottom left

	residual.set( w/2-1, h/2-1, Cplx.Float.one());    // zero order

	// TODO: check if this is what matlabs conv2 would do
	Transforms.swapQuadrant( residual );
	residual.fft2d(true);

	spec.times(residual);
	spec.fft2d(false);

    }

    /** Helper function, returns the magnitude of an input vector */
    Vec2d.Real magnitude( Vec2d.Cplx in ) {
	Vec2d.Real ret = Vec2d.createReal(in);
	ret.copyMagnitude(in);
	return ret;
    }



    /** ImageJ plugin run method */
    @Override
    public void run(String arg) {

	// calculate gratings 3.2 .. 3.8 pxl size,
	// allowing for 3 equi-distant phase shifts
	Tool.tell("-- Compute all candidates --");
	List<Grating> all = calcGrat(2.477,2.6, 3 );

	// collect pairs of 3 directions
	Tool.tell("-- Compute direction combinations --");
	List<Grating []> dirs = selectDirs(all, 3, 5*Math.PI/180.);
	Tool.tell(" Number of possible combinations: "+dirs.size());

	/*
	for ( Grating [] i : dirs ) {
	    System.out.println("---");
	    for ( Grating j : i )
		System.out.println( j.toString());
	} */

	// check if unwanted orders can be blocked by masking
	DisplayWrapper imgSpatial = new DisplayWrapper(fsPxl, fsPxl, "Spatial");
	DisplayWrapper imgFourier = new DisplayWrapper(fsPxl, fsPxl, "Fourier");
	
	Vec2d.Real gaussProfile = createGaussIllum( fsPxl/2.2, fsPxl);
	
	for (int i=0; i<Math.min( dirs.size(), 10); i++) {
	    fourierCheck( dirs.get(i), gaussProfile, 0.03, 20, null, imgFourier); 
	}



	/*
	
	Vec2d.Real sum = Vec2d.createReal( fsPxl, fsPxl );
	ImageDisplay.Marker [] maskRings = new ImageDisplay.Marker[6];
	
	for (int i=0; i<3; i++) {
	    
	    Vec2d.Real testPttr = Vec2d.createReal(512,512);
	    dirs.get(0)[i].writeToVector( testPttr );
	
	    imgSpatial.addImage(testPttr, "i: "+i);


	    Vec2d.Cplx fft = Vec2d.createCplx( testPttr );
	    fft.copy( testPttr );
	    multGauss( fft, 210);
	    //convSLMstructure( fft );
	    fft.fft2d(false);
		
	    fft.scal(1.f/(fsPxl*fsPxl));

	    ImageVector fftPttr = ImageVector.create(512,512);
	    fftPttr.copyMagnitude( fft );

	    Transforms.swapQuadrant( fftPttr );
	    sum.add(fftPttr );
	    
	    // compute the mask position
	    double [] mPos =  dirs.get(0)[i].peakPos(512);
	    ImageDisplay.Marker mp1 = new ImageDisplay.Marker( 
		256+mPos[0], 256+mPos[1], 10, 10, true);
	    ImageDisplay.Marker mm1 = new ImageDisplay.Marker( 
		256-mPos[0], 256-mPos[1], 10, 10, true);

	    maskRings[2*i] =  mp1;
	    maskRings[2*i+1] =  mm1;

	    imgFourier.addImage(fftPttr, "i: "+i, mp1, mm1);

	}
	
	imgFourier.addImage( sum, "Sum", maskRings );
	*/

	//imgSpatial.display();
	imgFourier.display();
	
    }


    
    /** main method */
    public static void main( String [] args ) {
	Grating_Search gs  = new Grating_Search();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	gs.run("");
    }


}




