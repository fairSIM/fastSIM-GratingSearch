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

/** Java port of fastSIM SLM grating search algorithm.
The original MATLAB code can be found here, please cite their
publication if you use this software to create SIM gratings:

https://github.com/nanoimaging/fastSIM_GratingSearchforSLM
*/
public class Grating_Search {

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


    // TODO: Implement the fourier plane unwanted orders checks here
    public void fourierCheck( List<Grating []> candidates ) {

    }









    
    /** main method */
    public static void main( String [] args ) {

	Grating_Search gs  = new Grating_Search();

	// calculate gratings 3.2 .. 3.8 pxl size,
	// allowing for 3 equi-distant phase shifts
	System.out.println("-- Compute all candidates --");
	List<Grating> all = gs.calcGrat(3.2,3.8, 3 );

	// collect pairs of 3 directions
	System.out.println("-- Compute direction combinations --");
	List<Grating []> dirs = gs.selectDirs(all, 3, 0.04);

	for ( Grating [] i : dirs ) {
	    System.out.println("---");
	    for ( Grating j : i )
		System.out.println( j.toString());
	}

	
    }






}




