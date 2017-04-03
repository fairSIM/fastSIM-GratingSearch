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
import ij.ImagePlus;
import javax.swing.JFileChooser;

import org.fairsim.linalg.Vec2d;
import org.fairsim.fiji.DisplayWrapper;
import org.fairsim.fiji.ImageVector;

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
	     IJ.getProperty("de.bio_photonics.gratingsearch.lastGratings")==null ||
	     IJ.getProperty("de.bio_photonics.gratingsearch.prefix")==null ||
	     IJ.getProperty("de.bio_photonics.gratingsearch.width")==null ||
	     IJ.getProperty("de.bio_photonics.gratingsearch.height")==null 
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
// 	gd.addNumericField("SLM width",1280,0);
// 	gd.addNumericField("SIM height",1024,0);
// 	gd.addStringField("file prefix","automatic", 30);
	gd.showDialog();
	if (gd.wasCanceled())
	    return;
    
	// get and check parameters
	int nr = (int)gd.getNextNumber()-1;
// 	int width   = (int)gd.getNextNumber();
// 	int height  = (int)gd.getNextNumber();
	String prefix = (String)IJ.getProperty("de.bio_photonics.gratingsearch.prefix");
	int width = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.width");
	int height = (Integer)IJ.getProperty("de.bio_photonics.gratingsearch.height");
//     OpenDialog directory = new ij.io.OpenDialog("Choose a directory"); 
// 	String path = directory.getDirectory();
	final String path;
	final JFileChooser fc = new JFileChooser();
    fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    fc.setAcceptAllFileFilterUsed(false);
    fc.setDialogTitle("choose directory");
    int returnVal = fc.showSaveDialog(null);
    if(returnVal == JFileChooser.APPROVE_OPTION) {
		path = fc.getSelectedFile().getAbsolutePath();
	}
    else path="";
	

	// generate pattern
	DisplayWrapper img = new DisplayWrapper(width, height, "Pattern");

	int ang=0;
	for ( Grating [] gr : gratList.get(nr) ) {
	    for (int pha =0; pha<nrPhases; pha++) {
			double phase = pha*Math.PI*2/nrPhases;
			ImageVector pttr = ImageVector.create(width,height);
			for ( Grating gri : gr ) {
				gri.writeToVector( pttr , phase );
				img.addImage( pttr, String.format("%s_wl%.0f_ang%d_pha:%d", prefix, gri.wavelength, ang, pha*360/nrPhases));
				if(returnVal == JFileChooser.APPROVE_OPTION) {
					ImagePlus tmpIP = new ImagePlus("test", pttr.img());
					IJ.save(tmpIP, String.format("%s/%s_wl%.0f_ang%d_pha%d.bmp", path, prefix, gri.wavelength, ang, pha*360/nrPhases));
				}
			}
	    }
	    ang++;
	}
	img.display();
    }
}
