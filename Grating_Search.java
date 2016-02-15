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


    // all defining features of the grating
    class GratingParam {
	int ax,ay,bx,by;// integer pxl size and inc. of grating
	double gratdir;	// orientation of grating in rad
	double gratper;	// grating period, in pxl
    }


    // main loop creating possible gratings
    public void calcGrat(
	final double gratPerMin, 
	final double gratPerMax ) {

	for (int ax = 0; ax < axmax; ax ++ ) {
	    System.out.println(String.format("DONE: %d / %d ",ax,axmax));

	    for( int ay = -aymax; ay < aymax; ay++ ) {

		double gratDir = Math.atan2(ay,ax);

		for ( int bx =  bxmin; bx < bxmax; bx++ )
		for ( int by = -bxmax; by < bymax; by++ ) {
		     
		    if (Math.abs(by)<bymin) continue;

		    // get the grating period
		    double gratPer = gratingPeriod(ax, ay, bx, by);
		    if ( gratPer < gratPerMin || gratPer > gratPerMax )
			continue;


		    System.out.println("gratPer: "+gratPer);



		}
	    }
	}





    }


    /** Calculate the grating period */
    double gratingPeriod( int ax, int ay, int bx, int by ) {
	double phi	= Math.atan2( ay, ax ) ;
	double theta	= Math.atan2( by, bx ) ;
	double res	= Math.hypot( bx, by ) * Math.abs( Math.sin( phi-theta ));
	return res;
    }

    
    /** main method */
    public static void main( String [] args ) {

	Grating_Search gs = new Grating_Search();
	gs.calcGrat(3,4);

    }



}




