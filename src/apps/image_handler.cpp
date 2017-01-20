/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <src/image.h>
#include <src/funcs.h>
#include <src/args.h>
#include <src/fftw.h>
#include <src/time.h>
#include <src/symmetries.h>

class image_handler_parameters
{
	public:
   	FileName fn_in, fn_out, fn_sel, fn_img, fn_sym, fn_sub, fn_mult, fn_div, fn_add, fn_subtract, fn_fsc, fn_adjust_power;
	int bin_avg, avg_first, avg_last, edge_x0, edge_xF, edge_y0, edge_yF, filter_edge_width, new_box;
	bool do_add_edge, do_flipXY, do_flipmXY, do_shiftCOM, do_stats;
	DOUBLE multiply_constant, divide_constant, add_constant, subtract_constant, threshold_above, threshold_below, angpix, new_angpix, lowpass, highpass, bfactor, shift_x, shift_y, shift_z;
   	int verb;
	// I/O Parser
	IOParser parser;

	Image<DOUBLE> Iout;
	Image<DOUBLE> Iop;
	MetaDataTable MD;

	// Image size
	int xdim, ydim, zdim;
	long int ndim;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
	    fn_in = parser.getOption("--i", "Input STAR file, image (.mrc) or movie/stack (.mrcs)");
	    fn_out = parser.getOption("--o", "Output name (overwrite input if empty; for STAR-input: insert this string before each image's extension)", "");

	    int cst_section = parser.addSection("image-by-constant operations");
	    multiply_constant = textToFloat(parser.getOption("--multiply_constant", "Multiply the image(s) pixel values by this constant", "1"));
	    divide_constant = textToFloat(parser.getOption("--divide_constant", "Divide the image(s) pixel values by this constant", "1"));
	    add_constant = textToFloat(parser.getOption("--add_constant", "Add this constant to the image(s) pixel values", "0."));
	    subtract_constant = textToFloat(parser.getOption("--subtract_constant", "Subtract this constant from the image(s) pixel values", "0."));
	    threshold_above = textToFloat(parser.getOption("--threshold_above", "Set all values higher than this value to this value", "999."));
	    threshold_below = textToFloat(parser.getOption("--threshold_below", "Set all values lower than this value to this value", "-999."));

	    int img_section = parser.addSection("image-by-image operations");
	    fn_mult = parser.getOption("--multiply", "Multiply input image(s) by the pixel values in this image", "");
	    fn_div = parser.getOption("--divide", "Divide input image(s) by the pixel values in this image", "");
	    fn_add = parser.getOption("--add", "Add the pixel values in this image to the input image(s) ", "");
	    fn_subtract = parser.getOption("--subtract", "Subtract the pixel values in this image to the input image(s) ", "");
	    fn_fsc = parser.getOption("--fsc", "Calculate FSC curve of the input image with this image ", "");
	    fn_adjust_power = parser.getOption("--adjust_power", "Adjust the power spectrum of the input image to be the same as this image ", "");

	    int four_section = parser.addSection("per-image operations");
	    do_stats = parser.checkOption("--stats", "Calculate per-image statistics?");
	    bfactor = textToFloat(parser.getOption("--bfactor", "Apply a B-factor (in A^2)", "0."));
	    lowpass = textToFloat(parser.getOption("--lowpass", "Low-pass filter frequency (in A)", "-1."));
	    highpass = textToFloat(parser.getOption("--highpass", "High-pass filter frequency (in A)", "-1."));
	    angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in A)", "1."));
	    new_angpix = textToFloat(parser.getOption("--rescale_angpix", "Scale input image(s) to this new pixel size (in A)", "-1."));
	    new_box = textToInteger(parser.getOption("--new_box", "Resize the image(s) to this new box size (in pixel) ", "-1"));
	    filter_edge_width = textToInteger(parser.getOption("--filter_edge_width", "Width of the raised cosine on the low/high-pass filter edge (in resolution shells)", "2"));
	    do_shiftCOM = parser.checkOption("--shift_com", "Shift image(s) to their center-of-mass (only on positive pixel values)");
	    shift_x = textToFloat(parser.getOption("--shift_x", "Shift images this many pixels in the X-direction", "0."));
	    shift_y = textToFloat(parser.getOption("--shift_y", "Shift images this many pixels in the Y-direction", "0."));
	    shift_z = textToFloat(parser.getOption("--shift_z", "Shift images this many pixels in the Z-direction", "0."));

	    int three_d_section = parser.addSection("3D operations");
	    fn_sym = parser.getOption("--sym", "Symmetrise 3D map with this point group (e.g. D6)", "");

	    int preprocess_section = parser.addSection("2D-micrograph (or movie) operations");
	    do_flipXY = parser.checkOption("--flipXY", "Flip the image(s) in the XY direction?");
	    do_flipmXY = parser.checkOption("--flipmXY", "Flip the image(s) image(s)in the -XY direction?");
	    do_add_edge = parser.checkOption("--add_edge", "Add a barcode-like edge to the micrograph/movie frames?");
	    edge_x0 = textToInteger(parser.getOption("--edge_x0", "Pixel column to be used for the left edge", "0"));
	    edge_y0 = textToInteger(parser.getOption("--edge_y0", "Pixel row to be used for the top edge", "0"));
	    edge_xF = textToInteger(parser.getOption("--edge_xF", "Pixel column to be used for the right edge", "4095"));
	    edge_yF = textToInteger(parser.getOption("--edge_yF", "Pixel row to be used for the bottom edge", "4095"));

	    int avg_section = parser.addSection("Movie-frame averaging options");
       	bin_avg = textToInteger(parser.getOption("--avg_bin", "Width (in frames) for binning average, i.e. of every so-many frames", "-1"));
    	avg_first = textToInteger(parser.getOption("--avg_first", "First frame to include in averaging", "-1"));
    	avg_last = textToInteger(parser.getOption("--avg_last", "Last frame to include in averaging", "-1"));

    	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

    	verb = (do_stats || fn_fsc !="") ? 0 : 1;

	}



	void perImageOperations(Image<DOUBLE> &Iin, FileName &my_fn_out)
	{

		Image<DOUBLE> Iout;
		Iout().resize(Iin());

		if (do_add_edge)
		{
			// Treat X-boundaries
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				if (j < edge_x0)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), i, edge_x0);
				else if (j > edge_xF)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), i, edge_xF);
			}
			// Treat Y-boundaries
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				if (i < edge_y0)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), edge_y0, j);
				else if (i > edge_yF)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), edge_yF, j);
			}
		}

		// Flipping: this needs to be done from Iin to Iout (i.e. can't be done on-line on Iout only!)
		if (do_flipXY)
		{
			// Flip X/Y
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				DIRECT_A2D_ELEM(Iout(), i, j) = DIRECT_A2D_ELEM(Iin(), j, i);

			}
		}
		else if (do_flipmXY)
		{
			// Flip mX/Y
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				DIRECT_A2D_ELEM(Iout(), i, j) = DIRECT_A2D_ELEM(Iin(), XSIZE(Iin()) - 1 - j, YSIZE(Iin()) - 1 - i);
			}
		}
		else
		{
			Iout = Iin;
		}

		// From here on also 3D options
		if (fabs(multiply_constant - 1.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) *= multiply_constant;
			}
		}
		else if (fabs(divide_constant - 1.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) /= divide_constant;
			}
		}
		else if (fabs(add_constant) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) += add_constant;
			}
		}
		else if (fabs(subtract_constant) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) -= subtract_constant;
			}
		}
		else if (fn_mult != "")
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) *= DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_div != "")
		{
                    bool is_first = true;
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
                            if (ABS(DIRECT_A3D_ELEM(Iop(), k, i, j)) < 1e-10)
                            {
				if (is_first)
                                {
                                    std::cout << "Warning: ignore very small pixel values in divide image..." << std::endl;
                                    is_first = false;
                                }
                                DIRECT_A3D_ELEM(Iout(), k, i, j) = 0.;
                            }
                            else
				DIRECT_A3D_ELEM(Iout(), k, i, j) /= DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_add != "")
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) += DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_subtract != "")
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) -= DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_fsc != "")
		{
			/// TODO
			MultidimArray<DOUBLE> fsc;
			MetaDataTable MDfsc;
			getFSC(Iout(), Iop(), fsc);
			MDfsc.setName("fsc");
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc)
			{
				MDfsc.addObject();
				DOUBLE res = (i > 0) ? (XSIZE(Iout()) * angpix / (DOUBLE)i) : 999.;
				MDfsc.setValue(EMDL_SPECTRAL_IDX, (int)i);
				MDfsc.setValue(EMDL_RESOLUTION, 1./res);
				MDfsc.setValue(EMDL_RESOLUTION_ANGSTROM, res);
				MDfsc.setValue(EMDL_POSTPROCESS_FSC_GENERAL, DIRECT_A1D_ELEM(fsc, i) );
			}
			MDfsc.write(std::cout);
		}
		else if (fn_adjust_power != "")
		{
			MultidimArray<DOUBLE> spectrum;
			getSpectrum(Iop(), spectrum, AMPLITUDE_SPECTRUM);
			adaptSpectrum(Iin(), Iout(), spectrum, AMPLITUDE_SPECTRUM);
		}

		if (fabs(bfactor) > 0.)
			applyBFactorToMap(Iout(), bfactor, angpix);

		if (lowpass > 0.)
			lowPassFilterMap(Iout(), lowpass, angpix, filter_edge_width);

		if (highpass > 0.)
			highPassFilterMap(Iout(), highpass, angpix, filter_edge_width);

		// Shifting
		if (do_shiftCOM)
			selfTranslateCenterOfMassToCenter(Iout(), DONT_WRAP);
		else if (fabs(shift_x) > 0. || fabs(shift_y) > 0. || fabs(shift_z) > 0.)
		{
			Matrix1D<DOUBLE> shift(2);
			XX(shift) = shift_x;
			YY(shift) = shift_y;
			if (zdim > 1)
			{
				shift.resize(3);
				ZZ(shift) = shift_z;
			}
			selfTranslate(Iout(), shift, DONT_WRAP);
		}

		// Re-scale
		if (new_angpix > 0.)
		{
			int oldsize = XSIZE(Iout());
			int newsize = ROUND(oldsize * (angpix / new_angpix));
			resizeMap(Iout(), newsize);
			// Also reset the sampling rate in the header
			Iout.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, new_angpix);
			Iout.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, new_angpix);
			if (Iout().getDim() == 3)
				Iout.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, new_angpix);
		}
		// Re-window
		if (new_box > 0)
		{
			Iout().setXmippOrigin();
			if (Iout().getDim() == 2)
			{
				Iout().window(FIRST_XMIPP_INDEX(new_box), FIRST_XMIPP_INDEX(new_box),
						   LAST_XMIPP_INDEX(new_box),  LAST_XMIPP_INDEX(new_box));
			}
			else if (Iout().getDim() == 3)
			{
				Iout().window(FIRST_XMIPP_INDEX(new_box), FIRST_XMIPP_INDEX(new_box), FIRST_XMIPP_INDEX(new_box),
						   LAST_XMIPP_INDEX(new_box),  LAST_XMIPP_INDEX(new_box),  LAST_XMIPP_INDEX(new_box));
			}
		}

		if (fn_sym != "")
			symmetriseMap(Iout(), fn_sym);

		// Thresholding (can be done after any other operation)
		if (fabs(threshold_above - 999.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				if (DIRECT_A3D_ELEM(Iout(), k, i, j) > threshold_above)
					DIRECT_A3D_ELEM(Iout(), k, i, j) = threshold_above;
			}
		}
		if (fabs(threshold_below + 999.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				if (DIRECT_A3D_ELEM(Iout(), k, i, j) < threshold_below)
					DIRECT_A3D_ELEM(Iout(), k, i, j) = threshold_below;
			}
		}

		// Write out the result
		// Check whether fn_out has an "@": if so REPLACE the corresponding frame in the output stack!
		long int n;
		FileName fn_tmp;
		my_fn_out.decompose(n, fn_tmp);
		n--;
		if (n >= 0) // This is a stack...
		{

			// If an output name was not specified: just replace the input image (in an existing stack)
			if (fn_out == "")
			{
				Iout.write(fn_tmp, n, true, WRITE_REPLACE); // replace an image in an existing stack
			}
			else
			{
				// The following assumes the images in the stack come ordered...
				if (n == 0)
					Iout.write(fn_tmp, n, true, WRITE_OVERWRITE); // make a new stack
				else
					Iout.write(fn_tmp, n, true, WRITE_APPEND);
			}
		}
		else
			Iout.write(my_fn_out);

	}

	void run()
	{

		// By default: write single output images

		// Get a MetaDataTable
		if (fn_in.getExtension() == "star")
		{
			MD.read(fn_in);
		}
		else if (fn_in.getExtension() == "mrcs" && !fn_in.contains("@"))
		{
			if (bin_avg > 0 || (avg_first >= 0 && avg_last >= 0))
			{
				MD.addObject();
				MD.setValue(EMDL_IMAGE_NAME, fn_in);
			}
			else
			{
				// Read the header to get the number of images inside the stack and generate that many lines in the MD
				Image<DOUBLE> tmp;
				FileName fn_tmp;
				tmp.read(fn_in, false); //false means do not read image now, only header
				for (int i = 1; i <= NSIZE(tmp()); i++)
				{
					MD.addObject();
					fn_tmp.compose(i, fn_in);
					MD.setValue(EMDL_IMAGE_NAME, fn_tmp);
				}
			}
		}
		else
		{
			// Just individual image input
			MD.addObject();
			MD.setValue(EMDL_IMAGE_NAME, fn_in);
		}

		int i_img = 0;
		time_config();
   		if (verb > 0)
   			init_progress_bar(MD.numberOfObjects());

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			FileName fn_img;
			MD.getValue(EMDL_IMAGE_NAME, fn_img);

			Image<DOUBLE> Iin;
			// Initialise for the first image
			if (i_img == 0)
			{
				Image<DOUBLE> Ihead;
				Ihead.read(fn_img, false);
				Ihead.getDimensions(xdim, ydim, zdim, ndim);

				if (zdim > 1 && (do_add_edge || do_flipXY || do_flipmXY))
					REPORT_ERROR("ERROR: you cannot perform 2D operations like --add_edge, --flipXY or --flipmXY on 3D maps. If you intended to operate on a movie, use .mrcs extensions for stacks!");

				if (zdim > 1 && (bin_avg > 0 || (avg_first >= 0 && avg_last >= 0)))
					REPORT_ERROR("ERROR: you cannot perform movie-averaging operations on 3D maps. If you intended to operate on a movie, use .mrcs extensions for stacks!");

				if (fn_mult != "")
					Iop.read(fn_mult);
				else if (fn_div != "")
					Iop.read(fn_div);
				else if (fn_add != "")
					Iop.read(fn_add);
				else if (fn_subtract != "")
					Iop.read(fn_subtract);
				else if (fn_fsc != "")
					Iop.read(fn_fsc);
				else if (fn_adjust_power != "")
					Iop.read(fn_adjust_power);

				if (fn_mult != "" || fn_div != "" || fn_add != "" || fn_subtract != "" || fn_fsc != "" || fn_adjust_power != "")
					if (XSIZE(Iop()) != xdim || YSIZE(Iop()) != ydim || ZSIZE(Iop()) != zdim)
						REPORT_ERROR("Error: operate-image is not of the correct size");

			}


			if (do_stats) // only write statistics to screen
			{
				Iin.read(fn_img);
				DOUBLE avg, stddev, minval, maxval;
				Iin().computeStats(avg, stddev, minval, maxval);
				std::cout << fn_img << " : (x,y,z,n)= " << XSIZE(Iin()) << " x "<< YSIZE(Iin()) << " x "<< ZSIZE(Iin()) << " x "<< NSIZE(Iin()) << " ; avg= " << avg << " stddev= " << stddev << " minval= " <<minval << " maxval= " << maxval << std::endl;
			}
			else if (bin_avg > 0 || (avg_first >= 0 && avg_last >= 0))
			{
				// movie-frame averaging operations
				int avgndim = 1;
				if (bin_avg > 0)
				{
					avgndim = ndim / bin_avg;
				}
				Image<DOUBLE> Iavg(xdim, ydim, zdim, avgndim);

				if (ndim == 1)
					REPORT_ERROR("ERROR: you are trying to perform movie-averaging options on a single image/volume");

				FileName fn_ext = fn_out.getExtension();
				if (NSIZE(Iavg()) > 1 && ( fn_ext.contains("mrc") && !fn_ext.contains("mrcs") ) )
					REPORT_ERROR("ERROR: trying to write a stack into an MRC image. Use .mrcs extensions for stacks!");

				for (long int nn = 0; nn < ndim; nn++)
				{
					Iin.read(fn_img, true, nn);
					if (bin_avg > 0)
					{
						int myframe = nn / bin_avg;
						if (myframe < avgndim)
						{
							FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
							{
								DIRECT_NZYX_ELEM(Iavg(),myframe,0,i,j) += DIRECT_A2D_ELEM(Iin(), i, j); // just store sum
							}
						}
					}
					else if (avg_first >= 0 && avg_last >= 0 && nn+1 >= avg_first && nn+1 <= avg_last) // add one to start counting at 1
					{
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iin())
						{
							DIRECT_MULTIDIM_ELEM(Iavg(), n) += DIRECT_MULTIDIM_ELEM(Iin(), n); // just store sum
						}
					}
				}
				Iavg.write(fn_out);
			}
			else
			{
				Iin.read(fn_img);
				FileName my_fn_out;
				if (fn_out == "")
					my_fn_out = fn_img;
				else
				{
					if (fn_in.getExtension() == "star" || (fn_in.getExtension() == "mrcs" && !fn_in.contains("@")))
					{
						my_fn_out = fn_img.insertBeforeExtension("_" + fn_out);
					}
					else
					{
						my_fn_out = fn_out;
					}
				}
				perImageOperations(Iin, my_fn_out);
			}

			i_img++;
			if (verb > 0)
				progress_bar(i_img);
		}

		if (verb > 0)
			progress_bar(MD.numberOfObjects());

	}


};


int main(int argc, char *argv[])
{
	image_handler_parameters prm;

	try
    {

		prm.read(argc, argv);

		prm.run();

    }
    catch (RelionError XE)
    {
        prm.usage();
        std::cout << XE;
        exit(1);
    }
    return 0;
}



