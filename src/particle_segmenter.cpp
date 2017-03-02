/***************************************************************************
 *
 * Author: "Qiang Zhou"
 * School of Life Sciences, Tsinghua University, Beijing, China
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
#include "src/particle_segmenter.h"
#include "src/symmetries.h"
//#define DEBUG

void ParticleSegmenter::read(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Input STAR file ");
	fn_ref = parser.getOption("--ref", "STAR file with the reference names, or an MRC stack file.");
	fn_out = parser.getOption("--o", "Output rootname", "seg");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms"));
	particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)"));
	lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter in Angstroms for the references (prevent Einstein-from-noise!)","-1"));
	lowpass_particle = textToFloat(parser.getOption("--lowpass_particle", "Lowpass filter in Angstroms for the particles","-1"));
	do_invert = parser.checkOption("--invert", "Density in particles is inverted w.r.t. density in template");
	// Modified by ZhouQ
	int refsub_section = parser.addSection("Subtract reference options");
	do_outputrefsubptcls = parser.checkOption("--outputrefsubptcls", "Output ref-subtracted particles?");
	fn_refsubptcls = parser.getOption("--refsubptcls", "Root name of output ref-subtracted particles","refsub");
	//do_fitleastsquare = parser.checkOption("--fitleastsq", "Fit particles with references by least square method?");
	//do_fitscale = parser.checkOption("--fitscale", "Fit particles scales?");
       	do_rewindowrefsub = parser.checkOption("--rewindowrefsub", "Output rewindowed ref-subtracted particles?");
	do_rewindowraw = parser.checkOption("--rewindowrawptcls", "Output rewindowed raw particles?");
	fn_winptcls = parser.getOption("--rewinedptcls", "Root name of output rewindowed raw particles","rewin");
        rewindow_size = textToInteger(parser.getOption("--rewindow_size", "Size of rewindowed particles","-1"));
	do_normalise = parser.checkOption("--norm", "Normalise the background to average zero and stddev one?");
	do_ramp = !parser.checkOption("--no_ramp", "Just subtract the background mean in the normalisation, instead of subtracting a fitted ramping background. ");
        white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
        black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));
	new_particle_diameter = textToFloat(parser.getOption("--new_particle_diameter", "Diameter of the re-windowed particles in Angstroms","-1."));
	do_shiftinfouriertransform = parser.checkOption("--shift_in_fourier", "Shift particles in Fourier transform?");
        t_x = textToFloat(parser.getOption("--target_x", "Coordinate X of target center in pixel. The origin is the center of reference. ","0."));
        t_y = textToFloat(parser.getOption("--target_y", "Coordinate Y of target center in pixel. The origin is the center of reference. ","0."));
        t_z = textToFloat(parser.getOption("--target_z", "Coordinate Z of target center in pixel. The origin is the center of reference. ","0."));
	t_rot = textToFloat(parser.getOption("--target_rot", "AngleRot of current target relative to reference target.","0."));
	t_tilt = textToFloat(parser.getOption("--target_tilt", "AngleTilt of current target relative to reference target.","0."));
	t_psi = textToFloat(parser.getOption("--target_psi", "AnglePsi of current target relative to reference target.","0."));
	fn_symin = parser.getOption("--symin", "The present symmetry from which expand the euler angles distribution to c1 symmetry","c1");
	do_randomsym = parser.checkOption("--randomsym", "Expand (Randomize) current symmetry to C1 symmetry. Overrided by index_sym.");
	index_sym = textToInteger(parser.getOption("--index_sym", "Given index of symmetry to be transformed to. Starts from 0. Nothing to do if less than 0. Override randomsym","-1"));

	do_ctf = parser.checkOption("--ctf", "Perform CTF correction on the references?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
	do_onlymakestar = parser.checkOption("--onlymakestar", "Only output star file with updated location information of segmented particles, do not output reference-subtracted or rewindowed raw particles. This option overwrites do_rewindowraw and do_outputrefsubptcls.");
	remain_rawcoord = parser.checkOption("--remain_rawcoord", "Remain coordinates of raw particles while do_onlymakestar.");
	updateOffset2SegPart = parser.checkOption("--updateOffset2SegPart", "Update offset to segmented particle while do_onlymakestar.");
	fn_onlystar = parser.getOption("--suffix_onlystar", "Root name of particles when only output star file.","");

	int expert_section = parser.addSection("Expert options");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	//if (parser.checkForErrors())
	//	REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void ParticleSegmenter::usage()
{
	// Modified by ZhouQ
	// for Update info
	std::cout << "Update 2016/05/02" << std::endl;
	parser.writeUsage(std::cerr);
}

void ParticleSegmenter::initialise()
{
	if (do_outputrefsubptcls && do_rewindowraw)
	{
		REPORT_ERROR("ParticleSegmenter::initialise ERROR: can not output reference-subtracted particles and rewindow particles at the same time! Provide either of them.");
	}

	// Read in input metadata file
	MDin.read(fn_in);

	if (fn_out != "")
		fn_out = fn_in.withoutExtension() + "_" + fn_out;
	else
		REPORT_ERROR("ParticleSegmenter::initialise ERROR: Please provide the root name of the output star file.");
	// If new particle diameter is not set, set it to original diameter
	if ( new_particle_diameter < 0 )
	{
		new_particle_diameter = particle_diameter;
	}
	if (fn_symin != "c1" && fn_symin != "C1")
	{
		int pgGroup, pgOrder;
		SymList SL;
		SL.isSymmetryGroup(fn_symin, pgGroup, pgOrder);
		SL.read_sym_file(fn_symin);
		Matrix2D<DOUBLE>  L(4, 4), R(4, 4);
		Matrix2D<DOUBLE>  Identity(3,3);
		Identity.initIdentity();
		R_repository.clear();
		L_repository.clear();
		R_repository.push_back(Identity);
		L_repository.push_back(Identity);
		for (int isym = 0; isym < SL.SymsNo(); isym++)
		{
			SL.get_matrices(isym, L, R);
			R.resize(3, 3);
			L.resize(3, 3);
			R_repository.push_back(R);
			L_repository.push_back(L);
		}

	}

	// Fill tabulated sine and cosine tables
	tab_sin.initialise(5000);
	tab_cos.initialise(5000);

	if (do_outputrefsubptcls && ! do_onlymakestar )
	{
		// Read in the references
		Mrefs.clear();
		if (fn_ref.isStarFile())
		{
			MetaDataTable MDref;
			if (fn_ref.contains("_model.star"))
			{
				MDref.read(fn_ref, "model_classes");
			}
			else
			{
				// Just read normal STAR file with user-provided references
				MDref.read(fn_ref);
			}

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDref)
			{
				// Get all reference images and their names
				Image<DOUBLE> Iref;
				FileName fn_img;
				if (!MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_img))
				{
					if (!MDref.getValue(EMDL_IMAGE_NAME, fn_img)) // perhaps the reference are called rlnImageName?
						REPORT_ERROR("ParticleSegmenter::initialise ERROR: either provide rlnReferenceImage or rlnImageName in the reference STAR file!");
				}
				Iref.read(fn_img);
				Iref().setXmippOrigin();
				Mrefs.push_back(Iref());
			}
		}
		else
		{
			Image<DOUBLE> Iref;
			Iref.read(fn_ref);
			Iref().setXmippOrigin();
			Mrefs.push_back(Iref());
		}
		particle_size = XSIZE(Mrefs[0]);
		// Invert references if necessary
		if (do_invert)
		{
			for (int iref = 0; iref < Mrefs.size(); iref++)
			{
				Mrefs[iref] *= -1.;
			}
		}

		// Calculate the size of the FTs given the low-pass filter
		if (lowpass < 0.)
		{
			lowpass = 2. * angpix;
			current_size = particle_size;
		}
		else
		{
			current_size = 2 * ROUND(particle_size * angpix / lowpass);
		}

		// Calculate (downsized) Fourier transforms of the references
		PPref.clear();
		MultidimArray<DOUBLE> dummy;
		Projector projector(particle_size);
		for (long iref = 0; iref < Mrefs.size(); iref++)
		{
			projector.computeFourierTransformMap(Mrefs[iref], dummy, current_size);
			PPref.push_back(projector);
		}
	}
#ifdef DEBUG
	std::cerr << "Finishing initialise" << std::endl;
#endif
}

void ParticleSegmenter::run()
{


	long int nr_parts = MDin.numberOfObjects();
	rewinOriginOffsets.resize(nr_parts, NR_REWINDOWOFFSETS);

	int barstep;
	if (verb > 0)
	{
		std::cout << "Segmenting for all input particles..." << std::endl;
		init_progress_bar(nr_parts);
		barstep = XMIPP_MAX(1, nr_parts / 60);
	}

	long int ipart = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
	{

		if (verb > 0 && ipart % barstep == 0)
			progress_bar(ipart);

		segmentOneParticle(ipart);

	    ipart++;
	}

	if (verb > 0)
		progress_bar(nr_parts);

	write();
}

void ParticleSegmenter::segmentOneParticle(long int ipart)
{
	FourierTransformer transformer;
	Image<DOUBLE> img;
	MultidimArray<DOUBLE> Mref_rot;
	MultidimArray<Complex > Fimg, Fref;

	// Get the image
	FileName fn_img;
	MDin.getValue(EMDL_IMAGE_NAME, fn_img, ipart);
	// Do Not Work
	//img.read(fn_img,! do_onlymakestar);
	//if (! do_onlymakestar )
	//{
	//	img.read(fn_img);
	//}
	//else
	//{
	//	img().resize(100,100);
	//	img().initZeros();
	//}
	img.read(fn_img);
	if ( do_outputrefsubptcls)
	{
		if (XSIZE(img()) != particle_size)
		{
			std::cerr << " fn_img= " << fn_img << " XSIZE(img())= " << XSIZE(img()) << " reference size= " << particle_size << std::endl;
			REPORT_ERROR("ParticleSegmenter::segmentOneParticle ERROR: References and images do not have the same size!");
		}
	}

	Matrix1D<DOUBLE> offset(2);
	MDin.getValue(EMDL_ORIENT_ORIGIN_X, XX(offset), ipart);
	MDin.getValue(EMDL_ORIENT_ORIGIN_Y, YY(offset), ipart);
	//offset.selfROUND();

	// Low-pass filter the image if necessary
	if (lowpass_particle > 0. && ! do_onlymakestar)
	{
		DOUBLE radius = XSIZE(img()) * angpix / lowpass;
		radius -= WIDTH_FMASK_EDGEB / 2.;
		DOUBLE radius_p = radius + WIDTH_FMASK_EDGEB;
		transformer.FourierTransform(img(), Fimg);
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fimg)
		{
			DOUBLE r = sqrt((DOUBLE)(kp*kp + ip*ip + jp*jp));
			if (r < radius)
				continue;
			else if (r > radius_p)
				DIRECT_A3D_ELEM(Fimg, k, i, j) = 0.;
			else
			{
				DIRECT_A3D_ELEM(Fimg, k, i, j) *= 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGEB);
			}
		}
		transformer.inverseFourierTransform(Fimg, img());
	}

	// Get the transformation parameters and which reference this is
	DOUBLE psi = 0.;
	DOUBLE rot = 0.;
	DOUBLE tilt = 0.;
	MDin.getValue(EMDL_ORIENT_ROT, rot, ipart);
	MDin.getValue(EMDL_ORIENT_TILT, tilt, ipart);
	MDin.getValue(EMDL_ORIENT_PSI, psi, ipart);
	// Get the rotational matrix
	Matrix2D<DOUBLE> A;
	Euler_angles2matrix(rot, tilt, psi, A);
	
	Matrix1D<DOUBLE> refsuboffset(2);
	// For rewindowing
	if ( do_rewindowrefsub || do_rewindowraw || do_onlymakestar )
	{
		Matrix1D<DOUBLE> p_ptcls;
		Matrix1D<DOUBLE> p_ref = vectorR3(t_x,t_y,t_z) ;
		p_ptcls = A * p_ref * -1. ;
		DIRECT_A2D_ELEM(rewinOriginOffsets, ipart, REWINDOWOFFSETX) = XX(p_ptcls);
		DIRECT_A2D_ELEM(rewinOriginOffsets, ipart, REWINDOWOFFSETY) = YY(p_ptcls);
		XX(refsuboffset) = floor(XX(offset) + DIRECT_A2D_ELEM(rewinOriginOffsets,ipart,REWINDOWOFFSETX));
		YY(refsuboffset) = floor(YY(offset) + DIRECT_A2D_ELEM(rewinOriginOffsets,ipart,REWINDOWOFFSETY));
	}

	if (do_rewindowraw && ! do_onlymakestar)
	{
		Image<DOUBLE> outrawimg;
		outrawimg = img;
		outrawimg().setXmippOrigin();
		if (do_shiftinfouriertransform)
		{
			transformer.FourierTransform(outrawimg(), Fimg);
			shiftImageInFourierTransform(Fimg, Fimg,tab_sin,tab_cos, XSIZE(outrawimg()), XX(refsuboffset),YY(refsuboffset));
			transformer.inverseFourierTransform(Fimg,outrawimg());
		}
		else
		{
			selfTranslate(outrawimg(), refsuboffset, true);
		}
		rewindow(outrawimg,rewindow_size);
		if (do_normalise) 
		{
			long int bg_rawradius = ROUND(new_particle_diameter / (2. * angpix));
			normalise(outrawimg,bg_rawradius,white_dust_stddev,black_dust_stddev,do_ramp);
		}
		FileName fn_outrawimg = fn_img.insertBeforeExtension("_" + fn_winptcls);
		outrawimg.write(fn_outrawimg,-1,true,2);
	}
	if ( do_outputrefsubptcls  && ! do_onlymakestar)
	{
		int iref;
		if (!MDin.getValue(EMDL_PARTICLE_CLASS, iref, ipart) && PPref.size() > 1)
			REPORT_ERROR("ParticleSegmenter::segment: cannot find class number in input STAR file");

		iref -= 1; // start counting at 0 instead of 1
		if (iref >= PPref.size())
		{
			MetaDataTable MDt;
			MDt.addObject(MDin.getObject(ipart));
			MDt.write(std::cerr);
			REPORT_ERROR("Too large class number for the given number of references: " + integerToString(iref));
		}
		// Get the reference image in the right orientation
		Fref.resize(current_size, current_size/2 + 1);
		// Euler angle of Align
		PPref[iref].get2DFourierTransform(Fref, A, IS_NOT_INV);
		if (do_ctf)
		{
			CTF ctf;
			MultidimArray<DOUBLE> Fctf;
			ctf.read(MDin, MDin, ipart);
			Fctf.resize(Fref);
			//ctf.getFftwImage(Fctf, XSIZE(img()), YSIZE(img()), angpix, false, false, intact_ctf_first_peak, true);
			ctf.getFftwImage(Fctf, current_size, current_size, angpix, false, false, intact_ctf_first_peak, true);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
			{
				DIRECT_MULTIDIM_ELEM(Fref, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
			}
		}
		if (lowpass > 0.)
		{
			// Set the Fourier transform back to its original size...
			MultidimArray<Complex > Frefp;
			windowFourierTransform(Fref, Frefp, XSIZE(img()));
			Fref = Frefp;
		}

		Mref_rot.resize(img());
		shiftImageInFourierTransform(Fref, Fref,tab_sin,tab_cos,XSIZE(Mref_rot), -XX(offset),-YY(offset));
		transformer.inverseFourierTransform(Fref, Mref_rot);
		CenterFFT(Mref_rot, true);
/*
	if ( do_fitscale)
	{
		// Modified by ZhouQ
		long int radiusp = ROUND(particle_diameter/(2. * angpix));
		MultidimArray<DOUBLE> iMsk;
		iMsk.initZeros(Mref_rot);
		FOR_ALL_ELEMENTS_IN_ARRAY2D(iMsk)
		{
			long int idx = ROUND(sqrt(i*i + j*j));
			if (idx < radiusp + 1)
			{
				DIRECT_A2D_ELEM(iMsk, i, j) = 1.;
			}
		}
		// Calculate the optimal scale between the image and the reference:
		DOUBLE sumxa = 0.;
		DOUBLE suma2 = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mref_rot)
		{
			if (DIRECT_MULTIDIM_ELEM(iMsk, n) > 0.5)
			{
				sumxa += DIRECT_MULTIDIM_ELEM(Mref_rot, n) * DIRECT_MULTIDIM_ELEM(img(), n);
				suma2 += DIRECT_MULTIDIM_ELEM(Mref_rot, n) * DIRECT_MULTIDIM_ELEM(Mref_rot, n);
			}
		}
		if (suma2 < 1e-10)
			REPORT_ERROR("ERROR: Empty reference encountered!!");
		DOUBLE scale = sumxa/suma2;
		Mref_rot *= scale;
	}

	if (do_outputrefsubptcls && do_fitleastsquare)
	{
		fitParticle(img,Mref_rot);
	}
*/
		// Calculate the difference between the particle and the oriented reference
		img() -= Mref_rot;
		img().setXmippOrigin();
		Image<DOUBLE> outrefsubimg;
		outrefsubimg = img;
		if (do_rewindowrefsub)
		{
			if (do_shiftinfouriertransform)
			{
				transformer.FourierTransform(outrefsubimg(), Fimg);
				shiftImageInFourierTransform(Fimg, Fimg,tab_sin,tab_cos, XSIZE(outrefsubimg()), XX(refsuboffset),YY(refsuboffset));
				transformer.inverseFourierTransform(Fimg,outrefsubimg());
			}
			else
			{
				selfTranslate(outrefsubimg(), refsuboffset, true);
			}
			rewindow(outrefsubimg,rewindow_size);
		}
		if (do_normalise) 
		{
			long int bg_refsubradius;
			if (do_rewindowrefsub)
				bg_refsubradius = ROUND(new_particle_diameter / (2. * angpix ));
			else
				bg_refsubradius = ROUND(particle_diameter / (2. * angpix ));
			normalise(outrefsubimg,bg_refsubradius,white_dust_stddev,black_dust_stddev,do_ramp);
		}
		FileName fn_outimg = fn_img.insertBeforeExtension("_" + fn_refsubptcls);
		outrefsubimg.write(fn_outimg,-1,true,2);
	}
}

/*
// Does not work for small particles
void ParticleSegmenter::fitParticle(Image<DOUBLE> &img, MultidimArray<DOUBLE> &Mref_rot)
{
	// set same Xmipp origin to both ref and img
	img().setXmippOrigin();
	Mref_rot.setXmippOrigin();
	long int radiusp = ROUND(particle_diameter/(2. * angpix));
	MultidimArray<DOUBLE> count, sum, sum_img;
	sum.initZeros(radiusp + 1);
	count.initZeros(radiusp + 1);
	sum_img.initZeros(radiusp + 1);
	FOR_ALL_ELEMENTS_IN_ARRAY2D(img())
	{
		long int idx = ROUND(sqrt(i*i + j*j));
		if (idx < radiusp + 1)
		{
			sum(idx) += DIRECT_A2D_ELEM(Mref_rot, i, j);
			sum_img(idx) += DIRECT_A2D_ELEM(img(), i, j);
			count(idx) += 1.;
		}
	}

	for (long int i = 0; i < radiusp + 1; i++)
	{
		if (count(i) > 0.)
		sum(i) /= count(i);
		sum_img(i) /= count(i);
	}
	fit_point2D      onepoint;
	std::vector<fit_point2D> points;
	for (long int i=0; i < radiusp + 1; i++)
	{
		//onepoint.w = 1.;
		onepoint.w = count(i);
		onepoint.x = sum(i);
		onepoint.y = sum_img(i);
		points.push_back(onepoint);
	}
	DOUBLE slope_x, intercept_x, ccf_x;
	fitStraightLine(points,slope_x,intercept_x,ccf_x);
	Mref_rot *= slope_x;
	Mref_rot += intercept_x;
}
*/

void ParticleSegmenter::write()
{

	FileName fn_img;
	FileName fn_tmp;

	// Modified by ZhouQ
	if ( do_outputrefsubptcls ||  do_rewindowraw || do_onlymakestar)
	{
		FileName fn_outimg;
		long int ipart = 0;
		DOUBLE xoffset,yoffset,coord_x,coord_y;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			if (do_rewindowrefsub || do_rewindowraw || (do_onlymakestar && ! remain_rawcoord))
			{
				MDin.getValue(EMDL_ORIENT_ORIGIN_X,xoffset,ipart);
				MDin.getValue(EMDL_ORIENT_ORIGIN_Y,yoffset,ipart);
				xoffset += DIRECT_A2D_ELEM(rewinOriginOffsets,ipart,REWINDOWOFFSETX);
				yoffset += DIRECT_A2D_ELEM(rewinOriginOffsets,ipart,REWINDOWOFFSETY);
				MDin.getValue(EMDL_IMAGE_COORD_X,coord_x,ipart);
				MDin.getValue(EMDL_IMAGE_COORD_Y,coord_y,ipart);
				coord_x -= floor(xoffset);
				coord_y -= floor(yoffset);
				MDin.setValue(EMDL_IMAGE_COORD_X,coord_x,ipart);
				MDin.setValue(EMDL_IMAGE_COORD_Y,coord_y,ipart);
				MDin.setValue(EMDL_ORIENT_ORIGIN_X,0.,ipart);
				MDin.setValue(EMDL_ORIENT_ORIGIN_Y,0.,ipart);
			}
			if (do_onlymakestar && remain_rawcoord && updateOffset2SegPart )
			{
				MDin.getValue(EMDL_ORIENT_ORIGIN_X,xoffset,ipart);
				MDin.getValue(EMDL_ORIENT_ORIGIN_Y,yoffset,ipart);
				xoffset += DIRECT_A2D_ELEM(rewinOriginOffsets,ipart,REWINDOWOFFSETX);
				yoffset += DIRECT_A2D_ELEM(rewinOriginOffsets,ipart,REWINDOWOFFSETY);
				MDin.setValue(EMDL_ORIENT_ORIGIN_X,xoffset,ipart);
				MDin.setValue(EMDL_ORIENT_ORIGIN_Y,yoffset,ipart);
			}
			if (fn_symin != "c1" && fn_symin != "C1")
			{
				DOUBLE psi, rot, tilt, psip,rotp,tiltp;
				MDin.getValue(EMDL_ORIENT_ROT, rot, ipart);
				MDin.getValue(EMDL_ORIENT_TILT, tilt, ipart);
				MDin.getValue(EMDL_ORIENT_PSI, psi, ipart);
				if (do_randomsym && index_sym < 0)
				{
					// Rondom choose one of symmetry-related Euler angles
					index_sym = (int) floor(rnd_unif() * R_repository.size());
					if ( index_sym == R_repository.size() )
					{
						index_sym = 0;
					}
				}
				if (index_sym >= 0)
				{
					Euler_apply_transf(L_repository[index_sym], R_repository[index_sym], rot, tilt, psi, rotp, tiltp, psip);
				}
				MDin.setValue(EMDL_ORIENT_ROT, rotp);
				MDin.setValue(EMDL_ORIENT_TILT, tiltp);
				MDin.setValue(EMDL_ORIENT_PSI, psip);
			}

			if (t_rot != 0. || t_tilt != 0. || t_psi != 0.)
			{
				DOUBLE psi, rot, tilt, psip,rotp,tiltp;
				MDin.getValue(EMDL_ORIENT_ROT, rot, ipart);
				MDin.getValue(EMDL_ORIENT_TILT, tilt, ipart);
				MDin.getValue(EMDL_ORIENT_PSI, psi, ipart);
				// set an aux Identity matrix
				Matrix2D<DOUBLE>  Identity(3,3);
				Identity.initIdentity();
				// Generate matrix from the Euler angles of target
				Matrix2D<DOUBLE> euler(3, 3), t_euler;
				Euler_angles2matrix(t_rot, t_tilt, t_psi, t_euler);
				// Right matrix
				Euler_apply_transf(Identity, t_euler, rot, tilt, psi, rotp, tiltp, psip);
				MDin.setValue(EMDL_ORIENT_ROT, rotp);
				MDin.setValue(EMDL_ORIENT_TILT, tiltp);
				MDin.setValue(EMDL_ORIENT_PSI, psip);
			}

			MDin.getValue(EMDL_IMAGE_NAME, fn_img,ipart);
			if (do_outputrefsubptcls)
			{
				fn_outimg = fn_img.insertBeforeExtension("_" + fn_refsubptcls);
				MDin.setValue(EMDL_IMAGE_NAME, fn_outimg);
			}
			if (do_rewindowraw)
			{
				fn_outimg = fn_img.insertBeforeExtension("_" + fn_winptcls);
				MDin.setValue(EMDL_IMAGE_NAME, fn_outimg);
			}
			if (do_onlymakestar)
			{
				if (fn_onlystar != "" )
				{
					fn_outimg = fn_img.insertBeforeExtension("_" + fn_onlystar);
					MDin.setValue(EMDL_IMAGE_NAME, fn_outimg);
				}
			}
			ipart++;
		}
		if (do_outputrefsubptcls)
		{
			fn_tmp = fn_out + "_" + fn_refsubptcls + ".star";
			MDin.write(fn_tmp);
		}
		if (do_rewindowraw)
		{
			fn_tmp = fn_out + "_" + fn_winptcls +".star";
			MDin.write(fn_tmp);
		}
		if (do_onlymakestar)
		{
			fn_tmp = fn_out + ".star";
			MDin.write(fn_tmp);
		}
	}

	if (verb>0)
		std::cout <<" Written out "<< fn_tmp << std::endl;
}



