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

import java.util.List;
import java.util.ArrayList;

import ij.IJ;
import ij.gui.GenericDialog;

import org.fairsim.linalg.Vec2d;
import org.fairsim.fiji.DisplayWrapper;

public class Generate_Grating implements ij.plugin.PlugIn {

    public void run(String arg) {

	if (arg=="clear") {
	    IJ.log("clearing stored gratings");
	    IJ.setProperty("de.bio_photonics.gratingsearch.phaseNumber",null);
	    IJ.setProperty("de.bio_photonics.gratingsearch.lastGratings",null);
	    return;
	}

	// check / retrieve parameters
	if ( IJ.getProperty("de.bio_photonics.gratingsearch.phaseNumber")==null ||
	     IJ.getProperty("de.bio_photonics.gratingsearch.lastGratings")==null 
	    ) {
	    IJ.log("No pattern information stored, search for pattern first!");
	    return;
	}

	// TODO: Check for cast errors?
	int nrPhases = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.phaseNumber");

	Object a= IJ.getProperty("de.bio_photonics.gratingsearch.lastGratings");
	List tmp = ( a instanceof List )?((List)a):(null);
	List<Grating [][]> gratList = new ArrayList<Grating [][]>();

	for ( Object b : tmp ) 
	    if ( b instanceof Grating[][] ) gratList.add((Grating [][])b);    
	

	if (gratList !=null && gratList.size()==0) {
	    IJ.log("No pattern in list, run search again");
	    return;
	}

	// show dialog
	GenericDialog gd = new GenericDialog("Pattern generation");
	gd.addNumericField(String.format("Pattern nr [1 - %d]",gratList.size()),1,0);
	gd.addNumericField("SLM width",1280,0);
	gd.addNumericField("SIM height",1024,0);

	gd.showDialog();
	if (gd.wasCanceled())
	    return;
    
	// get and check parameters
	int nr = (int)gd.getNextNumber()-1;
	int width = (int)gd.getNextNumber();
	int height= (int)gd.getNextNumber();

	if (nr<0 || nr>=gratList.size() || width<0 || height <0) {
	    IJ.log("Parameters out of range");
	    return;
	}

	// generate pattern
	DisplayWrapper img = new DisplayWrapper(width, height, "Pattern");

	int ang=0;
	for ( Grating [] gr : gratList.get(nr) ) {
	    for (int pha =0; pha<nrPhases; pha++) {
		double phase = pha*Math.PI*2/nrPhases;

		Vec2d.Real pttr = Vec2d.createReal(width,height);
		for ( Grating gri : gr ) {
		    gri.writeToVector( pttr , phase );
		    img.addImage( pttr, String.format("wl: %4.0f a: %d, pha: %d", gri.wavelength, ang, pha*360/nrPhases));
		}
	    }
	    ang++;
	}
	img.display();

    }

}
