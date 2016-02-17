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

import org.fairsim.linalg.Vec2d;
import org.fairsim.linalg.Cplx;
import org.fairsim.linalg.MTool;
import org.fairsim.utils.SimpleMT;

/** Structure hold all parameters of a grating. 
The original MATLAB code can be found here, please cite their
publication if you use this software to create SIM gratings:

https://github.com/nanoimaging/fastSIM_GratingSearchforSLM
*/
class Grating {

    final int ax,ay,bx,by;	// integer pxl size and inc. of grating
    final double gratDir;	// orientation of grating in rad
    final double gratPer;	// grating period, in pxl
    int shiftDir = -1;		// 0-horizontal, 1-vertical, -1 undef.

    /** create and calculate grating from pxl sizes */
    Grating(int iax, int iay, int ibx, int iby) {
	ax = iax; ay = iay; bx = ibx; by = iby;
	
	gratDir = Math.atan2(ay, ax);
	double theta	= Math.atan2( by, bx ) ;
	gratPer =	Math.hypot( bx, by ) * 
		    Math.abs( Math.sin( gratDir-theta ));
    }


    /** Test if a horizontal shift can provide n equi-distant
     *  phases. 
     *  @param nPhases Number of equi-distant phases to check */
    boolean testStepHorizontal( final int nPhases ) {

	// no shift for horizontal grating
	if (ay==0) return false;
	// standard shifts for vertical grating
	if (by==0) return ( bx % nPhases == 0 );

	int ax2 = ax/gcd(ax, ay);
	int ay2 = ay/gcd(ax, ay);

	return ( lcm(abs(ay2),abs(by))*(ax2/ay2-bx/by) % nPhases == 0 );

    }

    /** Test if a vertical shift can provide n equi-distant
     *  phases. 
     *  @param nPhases Number of equi-distant phases to check */
    boolean testStepVertical( final int nPhases ) {

	// no shift for horizontal grating
	if (ax==0) return false;
	// standard shifts for vertical grating
	if (bx==0) return ( by % nPhases == 0 );

	int ax2 = ax/gcd(ax, ay);
	int ay2 = ay/gcd(ax, ay);

	return ( lcm(abs(ax2),abs(bx))*(ay2/ax2-by/bx) % nPhases == 0 );
    }


    /** returns true if the pattern direction is within given range */
    boolean testDirection(double a1, double a2) {
	double d = wrapPi( gratDir );
	if (a1<=d && d<=a2) return true;
	return false;
    }


    @Override
    public String toString() {
	return String.format(
	    "ax %3d ay %3d bx %3d by %3d, d: %5.3f a %8.3f (%s)",
	    ax,ay,bx,by, gratPer, gratDir*180/Math.PI, 
	    (shiftDir==0)?("h"):("v"));
    }

    /** Write the (binary, -1, 1) pattern to a vector */
    public void writeToVector(final Vec2d.Real vec) {
	final double kx = (2*Math.PI / gratPer) * Math.sin(gratDir );
	final double ky = (2*Math.PI / gratPer) * Math.cos(gratDir );

	new SimpleMT.PFor(0, vec.vectorHeight()) {
	    public void at(int y) {
		for (int x=0; x<vec.vectorWidth(); x++) {
		    double val = MTool.fsin( Math.PI/2 + x*kx + y*ky + 1e-4);
		    vec.set(x,y,  (val>0)?(1):(-1));
		}
	    };
	};
    
    }

    /** Write the (binary, -1, 1) pattern to a vector */
    public void writeToVector(final Vec2d.Cplx vec) {
    
	final double kx = (2*Math.PI / gratPer) * Math.sin(gratDir );
	final double ky = (2*Math.PI / gratPer) * Math.cos(gratDir );

	new SimpleMT.PFor(0, vec.vectorHeight()) {
	    public void at(int y) {
		for (int x=0; x<vec.vectorWidth(); x++) {
		    double val = MTool.fsin( Math.PI/2 + x*kx + y*ky + 1e-4);
		    vec.set(x,y, new Cplx.Float( (val>0)?(1):(-1), 0 ));
		}
	    };
	};
    }
    
    /** Calculate the Fourier peak position caused by the grating,
     *  for a DFT of size n */
    public double [] peakPos(int n) {
	
	double len = n / gratPer;
	return new double [] { Math.sin(gratDir)*len, Math.cos(gratDir)*len };

    }


    /** Helper: greatest common denominator for int */
    static int gcd(int a, int b) { return b==0 ? a : gcd(b, a%b); }

    /** Helper: least common multiple */
    static int lcm(int a, int b) { 
	return (abs(a) / gcd(a,b))*abs(b); 
    }

    /** Helper: wrapper around abs */
    static int abs(int i) { return (i<0)?(-i):(i); };

    /** Helper: wrap direction into 0..pi */
    static double wrapPi(double x) {
	if (x>=Math.PI)	return wrapPi(x-Math.PI);
	if (x<0)	return wrapPi(x+Math.PI);
	return x;
    }


}


