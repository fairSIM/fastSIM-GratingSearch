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

import org.fairsim.utils.SimpleMT;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;

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
	    Tool.tell(String.format(
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
     *	@param	outputAlsoFailed Add also the gratings that failed unwanted order testing
     *	@param	spatial	    Receives spatial images of illuminated SLM, may be null
     *	@param  fourier	    Receives Fourier images of illuminated SLM, may be null
     *	@param name	    String to add to the displayed image caption
     * */
    public boolean fourierCheck( Grating [] candidates, Vec2d.Real illum, 
	double maxUnwated, int maskSize, boolean outputAlsoFailed,
	ImageDisplay spatial, ImageDisplay fourier, String name ) {

	Vec2d.Cplx     sumFreq = Vec2d.createCplx( fsPxl, fsPxl );
	Vec2d.Cplx [] gratFreq = Vec2d.createArrayCplx( candidates.length, fsPxl, fsPxl );

	// Compute the gratings Fourier space and sum them up
	for ( int i=0; i<candidates.length; i++ ) {
	    candidates[i].writeToVector( gratFreq[i] );	// get grating as vector 
	    gratFreq[i].times( illum );			// apply illumination vector
	    Transforms.fft2d( gratFreq[i], false);	// transform
	    Transforms.swapQuadrant(gratFreq[i]);	// quadrant-swap FFT result 
	    sumFreq.add( gratFreq[i] );			// add to full vector
	}

	// Compute wanted and unwanted contributions
	double []   wanted = new double[ candidates.length ];
	double [] unwanted = new double[ candidates.length ];
	boolean ok = true;
	String magResult = "";

	for ( int i=0; i<candidates.length; i++ ) {
	    
	    double [] mPos = candidates[i].peakPos(fsPxl);
	    
	    wanted[i] = 
		sumRegion( gratFreq[i],(int)(fsPxl/2+mPos[0]),(int)(fsPxl/2+mPos[1]), 8)+
		sumRegion( gratFreq[i],(int)(fsPxl/2-mPos[0]),(int)(fsPxl/2-mPos[1]), 8);

	    for (int j=0; j<candidates.length; j++) {
		if (i==j) continue;
		unwanted[i] += 
		   sumRegion( gratFreq[j],(int)(fsPxl/2+mPos[0]),(int)(fsPxl/2+mPos[1]), maskSize)+
		   sumRegion( gratFreq[j],(int)(fsPxl/2-mPos[0]),(int)(fsPxl/2-mPos[1]), maskSize);
	    }
	
	    //Tool.trace(String.format("%8.5f / %8.5f -> %8.5f", wanted[i], unwanted[i], 
	    //unwanted[i]/wanted[i]));
	    if ((unwanted[i]/wanted[i])>maxUnwated) ok=false;
	    magResult+=String.format("%1d: %5.4f ",i,unwanted[i]/wanted[i]);
	}

	// store spectrum (if display is set != null)
	if ( (fourier!=null) && (outputAlsoFailed || ok)) {
	    
	    // markers for positions
	    ImageDisplay.Marker [] maskRings = new ImageDisplay.Marker[candidates.length*4];
	    for ( int i=0; i<candidates.length; i++ ) {
		double [] mPos = candidates[i].peakPos(fsPxl);
		maskRings[i*4+0] = 
		    new ImageDisplay.Marker(fsPxl/2+mPos[0],fsPxl/2+mPos[1], 8, 8, false);
		maskRings[i*4+1] = 
		    new ImageDisplay.Marker(fsPxl/2-mPos[0],fsPxl/2-mPos[1], 8, 8, false);
		maskRings[i*4+2] = 
		    new ImageDisplay.Marker(fsPxl/2+mPos[0],fsPxl/2+mPos[1], maskSize, maskSize, false);
		maskRings[i*4+3] = 
		    new ImageDisplay.Marker(fsPxl/2-mPos[0],fsPxl/2-mPos[1], maskSize, maskSize, false);
	    }
	    fourier.addImage(magnitude( sumFreq ), 
		name+" "+magResult+((ok)?("ok"):("!!not OK")), maskRings);
	}

	// store pattern (if display is set != null)
	if ( (spatial != null) && (outputAlsoFailed || ok)) {
	    for (int i=0; i<candidates.length; i++) {
		Vec2d.Real img = Vec2d.createReal( fsPxl, fsPxl );
		candidates[i].writeToVector(img);
		img.times(illum);
		spatial.addImage( img, name+" ang: "+i+" "+candidates[i] );
	    }
	}
	
	return ok;
    }


    /** sum up the magnitude within a sub-region (mask) of a vector */
    double sumRegion( Vec2d.Cplx vec, int xp, int yp, int size ) {

	double sum=0;

	for (int y = Math.max(0, yp-size); y< Math.min( vec.vectorHeight()-1, yp+size); y++)
	for (int x = Math.max(0, xp-size); x< Math.min(  vec.vectorWidth()-1, xp+size); x++) {
	    sum += vec.get(x,y).abs();
	}

	return sum;
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

	Transforms.fft2d( spec, true);
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
	Transforms.fft2d( residual, true);

	spec.times(residual);
	Transforms.fft2d( spec, false);

    }

    /** Helper function, returns the magnitude of an input vector */
    Vec2d.Real magnitude( Vec2d.Cplx in ) {
	Vec2d.Real ret = Vec2d.createReal(in);
	ret.copyMagnitude(in);
	return ret;
    }



   
    /** Run the full calculation.
     *  @param gratMin	minial  grating period, in #pxl
     *  @param gratMax	maximal grating period, in #pxl
     *  @param phases	Number of phases
     *  @param nr_dir	Number of pattern orientations
     *  @param max_angle Maximal deviation from ideal (pi/n) angle distribution of pattern direction, in rad
     *  @param mask_size Size of holes in mask, in #pxl
     *  @param max_unwanted Maximal contribution of unwanted orders
     *  @param output_failed If to output also failed pattern
     *  */
    public void calculate(double gratMin, double gratMax, 
	int phases, int nr_dir, double max_angle, int mask_size, 
	double max_unwanted, boolean output_failed) {

	SimpleMT.useParallel(true);

	// calculate gratings 3.2 .. 3.8 pxl size,
	// allowing for 3 equi-distant phase shifts
	Tool.trace("-- Compute all candidates --");
	List<Grating> all = calcGrat(gratMin, gratMax, phases );
	Tool.trace(" Number matching grating constant and phases: "+all.size());

	// collect pairs of 3 directions
	Tool.trace("-- Compute direction combinations --");
	List<Grating []> dirs = selectDirs(all, nr_dir, max_angle);
	Tool.trace("   remaining number of possible combinations: "+dirs.size());

	/*
	for ( Grating [] i : dirs ) {
	    System.out.println("---");
	    for ( Grating j : i )
		System.out.println( j.toString());
	} */

	// check if unwanted orders can be blocked by masking
	Tool.trace("-- Compute wanted vs. unwanted orders --");
	DisplayWrapper imgSpatial = new DisplayWrapper(fsPxl, fsPxl, "Spatial");
	DisplayWrapper imgFourier = new DisplayWrapper(fsPxl, fsPxl, "Fourier");
	
	Vec2d.Real gaussProfile = createGaussIllum( fsPxl/2.2, fsPxl);
	int countOk=0, maxCheck = Math.min( dirs.size(), 400);

	for (int i=0; i< maxCheck ; i++) {
	    
	    if (i%10==0) Tool.tell(" FFT "+i+"/"+maxCheck);
	    
	    boolean ok = fourierCheck( dirs.get(i), gaussProfile, max_unwanted, mask_size, 
		output_failed, imgSpatial, imgFourier, "i:"+i); 
	    
	    if (ok) countOk++;
	}
	
	imgSpatial.display();
	imgFourier.display();
	Tool.trace("    Number of pattern after modulation check: "+ countOk);

    }

    // ---- Start methods ----

    /** ImageJ plugin run method */
    @Override
    public void run(String arg) {
	
	// redirect log output to FIJIs log
	Tool.setLogger( new Tool.Logger () {
	    @Override
	    public void writeTrace(String w) {
		    ij.IJ.log(w);
	    }
	    @Override
	    public void writeShortMessage(String w) {
	        ij.IJ.showStatus(w);
	    }

	});

	// generate simple GUI
	GenericDialog gd = new GenericDialog("SLM pattern search");

	gd.addMessage("Grating parameters");
	gd.addNumericField("grating_min_period", 2.81, 2);
	gd.addNumericField("grating_max_period", 2.85, 2);
	gd.addNumericField("#phases", 3, 0);
	gd.addMessage("Pattern direction");
	gd.addNumericField("#pattern_directions", 3, 0);
	gd.addNumericField("max_deviation_ideal_angle(deg)", 1.5, 1);
	gd.addMessage("Modulation");
	gd.addNumericField("max_unwanted_modulation",0.015,3);
	gd.addNumericField("mask_size(pxl)", 25, 0);
	gd.addCheckbox("Output_also_failed", false);
    
	gd.showDialog();
	if (gd.wasCanceled())
	    return;

	// get parameters
	double gratMin     = gd.getNextNumber();
	double gratMax     = gd.getNextNumber();
	int nrPhases  = (int)gd.getNextNumber();
	int nrDirs    = (int)gd.getNextNumber();
	double maxAngleDev = gd.getNextNumber();
	double maxUnwMod   = gd.getNextNumber();
	int maskSize  = (int)gd.getNextNumber();
	boolean outputFailed=gd.getNextBoolean();
	
	// run the actual calculation
	calculate( gratMin, gratMax, nrPhases, nrDirs, maxAngleDev, maskSize,
	    maxUnwMod, outputFailed);

	// output parameters at end of log
	Tool.tell(String.format("# Parameters: gMin %5.3f gMax %5.3f "+
	    " #Phases %1d  #Orientations %1d", gratMin, gratMax,
	    nrPhases, nrDirs));
	Tool.tell(String.format("# Search param: max angle deviation (deg): %6.0f",
	    maxAngleDev));
	Tool.tell(String.format("# Search param:   max unwanted modulation: %6.4f",
	    maxUnwMod));
	Tool.tell("# Search param: Mask size: "+maskSize);

    }
    

    
    /** main method */
    public static void main( String [] args ) {

	if (args.length!=1 && args.length!=4) {
	    System.out.println("Usage: gratPerMin gratPerMax #phases "+
		"#orientation");
	    System.out.println("or: GUI");
	    return;
	}
	
	Grating_Search gs  = new Grating_Search();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	if (args.length==1) {
	    gs.run("");
	    return;
	}

	double gratMin = Double.parseDouble( args[0] );
	double gratMax = Double.parseDouble( args[1] );
	int nrPhases   = Integer.parseInt( args[2] );
	int nrDirs     = Integer.parseInt( args[3] );

	gs.calculate(gratMin, gratMax, nrPhases, 
	    nrDirs, 3./180*Math.PI, 20, 0.02, false);
    }


}




