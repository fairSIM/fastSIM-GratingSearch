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
import ij.IJ;

/** Java port of fastSIM SLM grating search algorithm.
The original MATLAB code can be found here, please also cite their
publication if you use this software to create SIM gratings:

https://github.com/nanoimaging/fastSIM_GratingSearchforSLM
*/
public class Grating_Search implements ij.plugin.PlugIn {

    // search value space
    final int axmax = 30, aymax = 30;
    final int bxmin = 2,  bymin = 2;
    final int bxmax = 30, bymax = 30;

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
	final int phaseSteps, 
	final double wavelength ) {

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
		    Grating current = new Grating(ax,ay,bx,by, wavelength);
	   
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
	final double tolDir = 1e-8;
	
	for ( Grating b : candidates )
	    if ( Math.abs( b.gratDir  - a.gratDir ) < tolDir )
	    if ( Math.abs( b.gratPer  - a.gratPer ) < tolDir )
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
    public List<Grating [][]> selectDirs( 
	List< List<Grating> > candidates,
	final int nrWavelength, 
	final int nrDirs, 
	final double tolDirection,
	final double tolBetweenWavelenth,
	final int maxCandidates ) {

	double sector = Math.PI/nrDirs;

	int count=0;
	int maxCount = candidates.get(0).size();
    
	List<Grating [][]> ret = new ArrayList<Grating [][]>(maxCandidates);

	// outer loop, select only dirs in the first sector
	for ( Grating c : candidates.get(0) ) {
	    
	    count++;
	    if (count%100==0)
		ij.IJ.showStatus(String.format(" %d / %d tested, found: %d ",
		    count, maxCount, ret.size()));
	    
	    if (!c.testDirection(0, sector)) continue;

	    // create the array
	    Grating [][] col = new Grating[nrWavelength][nrDirs];
	    col[0][0] = c;


	    // see if the other directions can be found within limits
	    for (int i=1; i<nrDirs; i++) {
		double minDist=tolDirection;
		for ( Grating d : candidates.get(0) ) {
		    double dist =  Math.abs(
			    Grating.wrapPi(d.gratDir) - sector*i - 
			    Grating.wrapPi(c.gratDir) );
		    if ( dist < minDist ) {
			col[0][i]  = d;
			minDist = dist;
		    }
		}
	    }

	    // check if a complete set was found
	    boolean complete = true;
	    for ( Grating e : col[0] )
		if ( e == null ) continue;


	    // for all other wavelength, find the grating with minimal
	    // euclidian distance (corrected with wavelength to our candidate)
	    boolean allFound = true;
	    for (int  i=0; i<nrDirs; i++) {
		for (int ch=1; ch<nrWavelength; ch++) {
		    
		    double minDist = Double.MAX_VALUE;
		    for ( Grating d : candidates.get(ch) ) {
			double dist = Grating.scaledDistance( col[0][i], d );
			if (dist < minDist) {
			    minDist = dist;
			    col[ch][i] = d;
			}
		    }
		    // break the loop if no candidate was found within
		    // set distance limit
		    if (minDist > tolBetweenWavelenth)
			allFound = false;

		}
	    }
	    
	   
	    // only if all is good, add candidate to list
	    if (!allFound) continue;
	    ret.add( col );

	    if ( ret.size() > maxCandidates )
		break;
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
		candidates[i].writeToVector(img,0);
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

    // ugly trick to allow for array of generic typhes
    private interface ListOfGratings extends List<Grating []> {};


   
    /** Run the full calculation.
     *  @param wavelength list of wavelength to calculate for
     *  @param gratMin	for each wavelength, minial  grating period, in #pxl
     *  @param gratMax	for each wavelength, maximal grating period, in #pxl
     *  @param phases	Number of phases
     *  @param nr_dir	Number of pattern orientations
     *  @param max_angle Maximal deviation from ideal (pi/n) angle distribution of pattern direction, in rad
     *  @param mask_size Size of holes in mask, in #pxl
     *  @param max_unwanted Maximal contribution of unwanted orders
     *  @param output_failed If to output also failed pattern
     *  @param max_candidates How many candidates to calculate in first step
     *  @return List of matching gratings
     *  */
    public List<Grating [][]> calculate(
	double [] wavelength,
	double [] gratMin, double [] gratMax,
	int phases, int nr_dir, double max_angle, int mask_size, 
	double max_unwanted, double max_euclDist, 
	boolean output_failed, final int max_candidates) {

	List<Grating [][]> resultList = new ArrayList<Grating [][]>();
	
	SimpleMT.useParallel(true);

	Tool.trace("-- Compute all candidates --");

	// for each wavelength, create a list of candidates
	List< List<Grating> > all = new ArrayList< List<Grating> >();
	for (int ch=0; ch<wavelength.length; ch++) {
	    List< Grating > candidate = calcGrat( gratMin[ch], gratMax[ch], phases, wavelength[ch] );
	    all.add( candidate );
	    Tool.trace(String.format(" %5.0f nm : %d candidates", wavelength[ch], candidate.size()));
	}

	// check for candidates that overlap, so they fit through the same mask
	Tool.trace("-- Compute direction combinations --");
	List<Grating [][]> dirs = selectDirs(all, wavelength.length, nr_dir, 
	    max_angle, max_euclDist, max_candidates);
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
	int countOk=0;
	for (int i=0; i< dirs.size() ; i++) {
	    
	    if (i%10==0) Tool.tell(" FFT "+i+"/"+dirs.size());
	    
	    boolean ok = true;
	    
	    for (int ch=0; ch<wavelength.length; ch++) {
		ok &= fourierCheck( dirs.get(i)[ch], gaussProfile, max_unwanted, mask_size, 
		    output_failed, imgSpatial, imgFourier, "i:"+i); 
	    }
	    
	    if (ok) {
		countOk++;
		resultList.add( dirs.get(i) );
	    }
	}
	
	imgSpatial.display();
	imgFourier.display();
	Tool.trace("    Number of pattern after modulation check: "+ countOk);
	
	return resultList;
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

	gd.addMessage("System parameters");
	gd.addNumericField("proj. SLM pixel size (nm)", 70, 1);
	gd.addNumericField("resolution enhancement average  (at NA 1.45)", 1.75, 2);
	gd.addNumericField("resolution enhancement range +- (at NA 1.45)", 0.03, 2);

	gd.addMessage("Wavelength to analyse");
	gd.addCheckbox("include in calc", true);
	gd.addNumericField("wavelength 1: ", 488, 0);
	gd.addCheckbox("include in calc", true);
	gd.addNumericField("wavelength 2: ", 568, 0);
	gd.addCheckbox("include in calc", true);
	gd.addNumericField("wavelength 3: ", 647, 0);

	gd.addMessage("Pattern parameters");
	gd.addNumericField("#phases", 3, 0);
	gd.addNumericField("#pattern_directions", 3, 0);
	gd.addNumericField("max_deviation_ideal_angle(deg)", 1.5, 1);
	gd.addNumericField("max_eucl_dist(approx. pxl)", .05, 2);
	gd.addMessage("Modulation");
	gd.addNumericField("max_unwanted_modulation",0.015,3);
	gd.addNumericField("mask_size(pxl)", 15, 0);
	gd.addCheckbox("Output_also_failed", false);
	gd.addMessage("Cancel");
	gd.addNumericField("max_nr_candidates",500,0);

	gd.showDialog();
	if (gd.wasCanceled())
	    return;

	// ---- get parameters ----
	
	final double pxlSize	    = gd.getNextNumber();
	final double resImpAvr	    = gd.getNextNumber();
	final double resImpRange    = gd.getNextNumber();

	final double  [] wavelength	= new double[3];
	final boolean [] doWavelength	= new boolean[3];

	for (int i=0; i<3; i++) {
	    wavelength[i]   = gd.getNextNumber();
	    doWavelength[i] = gd.getNextBoolean();
	}
	
	final int nrPhases  = (int)gd.getNextNumber();
	final int nrDirs    = (int)gd.getNextNumber();
	final double maxAngleDev = gd.getNextNumber();
	final double maxEuclDist = gd.getNextNumber();
	
	final double maxUnwMod   = gd.getNextNumber();
	final int maskSize  = (int)gd.getNextNumber();
	final boolean outputFailed=gd.getNextBoolean();
	final int maxCandidates  = (int)gd.getNextNumber();
	
	
	// 1 - calculate the ranges for all from resolution enhancement
	//
	// This is pttr/2 = lambda / ( 2 * NA * (resImp-1) * pxlSize)
	// where
	//  pttr    : pattern spacing, in pxls
	//  lambda  : exitation wavelength
	//  NA	    : objective's NA
	//  resImp  : factor of resolution improvement, e.g. 1.75x
	//  pxlSize : projected pixel size 
	//
	final double objNA = 1.45;

	final double resImpMin = resImpAvr - resImpRange;
	final double resImpMax = resImpAvr + resImpRange;

	IJ.log(String.format("Searching pattern for res. improvement %7.4f -- %7.4f",
	    resImpMin, resImpMax ));

	double [] gratMin = new double[3];
	double [] gratMax = new double[3];
	

	for (int ch=0; ch<3; ch++) {
	    gratMax[ch] = 2 * wavelength[ch] / ( 2 * objNA * (resImpMin-1) * pxlSize );
	    gratMin[ch] = 2 * wavelength[ch] / ( 2 * objNA * (resImpMax-1) * pxlSize );

	    IJ.log(String.format("Search range %5.0f nm: %7.4f --- %7.4f pxl", wavelength[ch], 
		gratMin[ch], gratMax[ch]));

	}


	// run the actual calculation
	List<Grating [][]> res = calculate( wavelength, 
	    gratMin, gratMax, nrPhases, 
	    nrDirs, maxAngleDev, maskSize,
	    maxUnwMod, maxEuclDist, outputFailed, maxCandidates);

	/*
	// store the result, if any
	if (res.size()>0) {
	    ij.IJ.setProperty("de.bio_photonics.gratingsearch.phaseNumber", nrPhases);
	    ij.IJ.setProperty("de.bio_photonics.gratingsearch.lastGratings", res);
	}

	// output parameters at end of log
	Tool.tell(String.format("# Parameters: gMin %5.3f gMax %5.3f "+
	    " #Phases %1d  #Orientations %1d, max candidates: %d", gratMin, gratMax,
	    nrPhases, nrDirs, maxCandidates));
	Tool.tell(String.format("# Search param: max angle deviation (deg): %6.0f",
	    maxAngleDev));
	Tool.tell(String.format("# Search param:   max unwanted modulation: %6.4f",
	    maxUnwMod));
	Tool.tell("# Search param: Mask size: "+maskSize);
	*/

    }
    

    
    /** main method */
    public static void main( String [] args ) {

	
	Grating_Search gs  = new Grating_Search();
	new ij.ImageJ( ij.ImageJ.EMBEDDED);
	
	gs.run("");

	/*
	double gratMin = Double.parseDouble( args[0] );
	double gratMax = Double.parseDouble( args[1] );
	int nrPhases   = Integer.parseInt( args[2] );
	int nrDirs     = Integer.parseInt( args[3] );

	gs.calculate(gratMin, gratMax, nrPhases, 
	    nrDirs, 3./180*Math.PI, 20, 0.02, false, 400); */
    }


}




