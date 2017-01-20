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
#include "src/subtractingmic.h"

void SubtractingMicrograph::getMicrographsFromMetaData(MetaDataTable &MD, std::vector<FileName> &fn_allmics)
{
	fn_allmics.clear();
	FileName fn_mic;
	bool is_unique = true;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
		if ( do_rewindowParticles || (  do_rewindowMicrographs && exists(fn_mic)) || do_renumberRawImages )
		{
			is_unique = true;
			for (int imic = 0; imic < fn_allmics.size(); imic++)
			{
				if (fn_allmics[imic] == fn_mic)
				{
					is_unique = false;
					break;
				}
			}
			if (is_unique)
			{
				fn_allmics.push_back(fn_mic);
			}
		}
		else
		{
			if (do_rewindowMicrographs && verb > 0)
			{
				std::cout<<  "Warning: There is no micrograph:"<< fn_mic << std::endl;
			}
		}
	}
}


void SubtractingMicrograph::splitMetaDataAccordingtoMicrograph(MetaDataTable &MD, std::vector<MetaDataTable> &MD_mics, std::vector<FileName> &fn_mics, std::vector<CTF> &ctf_mics, std::vector<FileName> &fn_mymics)
{
	fn_mics.clear();
	ctf_mics.clear();
	MD_mics.clear();

	FileName fn_mic;
	bool is_unique = true;
	bool is_mine = false;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
		if ( fn_mymics.size() > 0 )
		{
			is_mine = false;
			for (int imic = 0; imic < fn_mymics.size(); imic++)
			{
				if (fn_mymics[imic] == fn_mic)
				{
					is_mine = true;
					break;
				}
			}
			if ( ! is_mine )
			{
				continue;
			}
		}
		if ( do_rewindowParticles || (  do_rewindowMicrographs && exists(fn_mic)) || do_renumberRawImages  ) 
		{
			is_unique = true;
			for (int imic = 0; imic < fn_mics.size(); imic++)
			{
				if (fn_mics[imic] == fn_mic)
				{
					is_unique = false;
					MD_mics[imic].addObject(MD.getObject());
					break;
				}
			}
			if (is_unique)
			{
				fn_mics.push_back(fn_mic);
				MetaDataTable MD_onemic;
				MD_onemic.addObject(MD.getObject());
				MD_mics.push_back(MD_onemic);
				// Store CTF for fn_mic
				CTF ctf;
				ctf.read(MD,MD,-1);
				ctf_mics.push_back(ctf);
			}
		}
		else
		{
			if (do_rewindowMicrographs && verb > 0)
			{
				std::cout<<  "Warning: There is no micrograph:"<< fn_mic << std::endl;
			}
		}
	}
}

void SubtractingMicrograph::read(int argc, char **argv,int rank)
{
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--input_star", "The input STAR file","");
	angpix  = textToFloat(parser.getOption("--angpix_micrograph", "Pixel size of RAW micrographs in Angstrom", "1.30654"));
	fn_coord = parser.getOption("--coord_star", "The STAR file providing coordinates of the particles rewindowing","");
	fn_sub = parser.getOption("--sub_star", "The STAR file to be subtracted","");
	fn_ref = parser.getOption("--ref", "The references file","");
	lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter (in Angstroms) for the input reference (default is none) prior to applying mask", "-1"));
	fn_mask = parser.getOption("--mask", "Mask applied to the input reference to select region to be projected and subtracted", "");
	binning  = textToFloat(parser.getOption("--binning", "Binning factor of references and input data relative to raw micrographs,e.g. 2. for bin2 ", "1."));
	do_invert_refs = parser.checkOption("--invert_refs", "Invert the contrast of references?");

	int renumraw_section = parser.addSection("Renumber raw image stack");
	do_renumberRawImages = parser.checkOption("--renumber_rawimages", "Renumber raw images?");
	fn_renumber = parser.getOption("--renumber_name", "The suffix for the renumbered star file","renumber");
	fn_renumber = fn_renumber.removeDirectories();

	int refsub_section = parser.addSection("Subtract projections or neighbor particles from micrographs");
	do_subtract = parser.checkOption("--subtract", "Subtract projection from micrographs or particles before re-extraction?");
	do_subtract_neighbor = parser.checkOption("--subtract_neighbor", "Subtract neighbor particles");
	neighbor_size = textToInteger(parser.getOption("--neighbor_size", "Pre-extracting box size before subtract neighbor particles  (in pixels, with same binning factor as reference, must be a even number.)", "1024"));
	do_output_refsub_micrographs = parser.checkOption("--output_refsub_micrographs", "Output refsub micrographs");
	fn_refsubmic = parser.getOption("--refsubmic", "The suffix for the refsub micrographs","refsub");
	fn_refsubmic = fn_refsubmic.removeDirectories();
	do_output_projection_micrographs = parser.checkOption("--output_projection_micrographs", "Output projection micrographs");
	fn_projectionmic = parser.getOption("--projectionmic", "The suffix for the projection micrographs","projection");
	fn_projectionmic = fn_projectionmic.removeDirectories();

	int ctf_section = parser.addSection("CTF section");
	do_ctf = parser.checkOption("--ctf", "Perform CTF correction on the projections?");
	do_apply_ctf_to_total = parser.checkOption("--apply_ctf_to_total_projection", "Perform CTF correction on the sum micrograph containing total projections?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");


	int reextract_section = parser.addSection("Re-extract particles from micrographs or particles");
	do_rewindowMicrographs = parser.checkOption("--rewindow_micrograph", "Rewindow micrographs?");
	do_invert_micrographs = parser.checkOption("--invert_micrographs", "Invert the contrast of micrographs?");
	do_rewindowParticles = parser.checkOption("--rewindow_particle", "Rewindow particles?");
	do_output_raw_particles = !parser.checkOption("--do_not_output_raw_particles", "Do not re-output raw particles matching with rewindowed particles.");
	do_output_sub_particles = !parser.checkOption("--do_not_output_sub_particles", "Do not re-output subtraction particles matching with rewindowed particles.");
	do_reextract = parser.checkOption("--reextract", "RE-extract particles when rewindowing micrographs?");
	fn_rewindow = parser.getOption("--extract_name", "The suffix for the rewindowed particles","rewindow");
	fn_rewindow = fn_rewindow.removeDirectories();
	extract_size = textToInteger(parser.getOption("--extract_size", "Size of the box to extract the particles in (in pixels, with same binning factor as reference)", "-1"));
	target_x = textToFloat(parser.getOption("--target_x", "Coordinate X of target relative to the centor of original map in pixel. ","0."));
	target_y = textToFloat(parser.getOption("--target_y", "Coordinate X of target relative to the centor of original map in pixel. ","0."));
	target_z = textToFloat(parser.getOption("--target_z", "Coordinate X of target relative to the centor of original map in pixel. ","0."));

	int perpart_section = parser.addSection("Particle operations");
	do_normalise = parser.checkOption("--norm", "Normalise the background to average zero and stddev one");
	do_ramp = !parser.checkOption("--no_ramp", "Just subtract the background mean in the normalisation, instead of subtracting a fitted ramping background. ");
	new_particle_diameter = textToFloat(parser.getOption("--new_particle_diameter", "Particle diameter of target part to be rewindowed in Angstrom", "100."));
	white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
	black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));
	// Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}
void SubtractingMicrograph::usage()
{
	std::cout << "Update 2016/11/30" << std::endl;
	parser.writeUsage(std::cerr);
}

void SubtractingMicrograph::initialise(int rank, int nodesize)
{
	if (fn_in == "")
	{
		REPORT_ERROR("ERROR: Please specify a input STAR file to be processed!");
	}
	if (do_subtract_neighbor && fn_coord != "")
	{
		REPORT_ERROR("ERROR: You can not specify coordinate STAR file during subtract neighbor particles.");
	}
	if (do_subtract_neighbor && (target_x != 0. || target_y != 0. || target_z != 0.))
	{
		REPORT_ERROR("ERROR: You can not specify target x/y/z during subtract neighbor particles.");
	}
	if (do_rewindowMicrographs && do_rewindowParticles)
	{
		REPORT_ERROR("ERROR: You can not specify --rewindow_micrograph and --rewindow_particle at the same time. Please choose one.");
	}
	// Find all unique micrographs in input STAR file
	MDin.read(fn_in);
	if (!MDin.containsLabel(EMDL_MICROGRAPH_NAME))
		REPORT_ERROR("ERROR: Input STAR file has no rlnMicrographName column!");
	long int my_first_mic, my_last_mic, my_nr_mics;
	fn_mymics.clear();
	if ( nodesize > 0 )
	{
		getMicrographsFromMetaData(MDin, fn_allmics);
		MDin.goToObject(0);
		divide_equally(fn_allmics.size(), nodesize, rank, my_first_mic, my_last_mic);
		my_nr_mics = my_last_mic - my_first_mic + 1;
		for (int imic = my_first_mic; imic <= my_last_mic; imic++)
		{
			fn_mymics.push_back(fn_allmics[imic]);
		}
	}
	splitMetaDataAccordingtoMicrograph(MDin,MD_mics,fn_mics,ctf_mics,fn_mymics);
	MDin.clear();

	// For coord STAR file
	if (fn_coord == "" )
	{
		if (verb > 0)
		{
			std::cout<<  "Warning: Using input star file: "<< fn_in << " as coordinate file." << std::endl;
		}
		fn_coord = fn_in;
	}
	// Find all unique micrographs in input STAR file
	MD_coord.read(fn_coord);
	if (!MD_coord.containsLabel(EMDL_MICROGRAPH_NAME))
		REPORT_ERROR("ERROR: Input coordinate STAR file has no rlnMicrographName column!");
	splitMetaDataAccordingtoMicrograph(MD_coord,MD_coord_mics,fn_coord_mics,ctf_coord_mics,fn_mymics);
	MD_coord.clear();
	
	// For subtract STAR file
	if (fn_sub == "" )
	{
		if (verb > 0)
		{
			std::cout<<  "Warning: Using input star file: "<< fn_in << " as subtract file." << std::endl;
		}
		fn_sub = fn_in;
	}
	// Find all unique micrographs in input STAR file
	MD_sub.read(fn_sub);
	if (!MD_sub.containsLabel(EMDL_MICROGRAPH_NAME))
		REPORT_ERROR("ERROR: Input subtract STAR file has no rlnMicrographName column!");
	splitMetaDataAccordingtoMicrograph(MD_sub,MD_sub_mics,fn_sub_mics,ctf_sub_mics,fn_mymics);
	MD_sub.clear();

	// Fill tabulated sine and cosine tables
	tab_sin.initialise(5000);
	tab_cos.initialise(5000);

	// Find all references
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
					REPORT_ERROR("subtract_micrograph::initialise ERROR: either provide rlnReferenceImage or rlnImageName in the reference STAR file!");
			}
			Iref.read(fn_img);
			Iref().setXmippOrigin();
			Mrefs.push_back(Iref());
		}
	}
	else
	{
		if (fn_ref != "")
		{
			Image<DOUBLE> Iref;
			Iref.read(fn_ref);
			Iref().setXmippOrigin();
			Mrefs.push_back(Iref());
		}
	}

	//Update ang pix of micrograph if references are scaled relative to micrographs.
	if (ABS(binning - 1.0) > 0.001 )
	{
		angpix *=  binning ;
	}

	if (lowpass > 0)
	{
		for(int iref=0; iref < Mrefs.size();iref++)
		{
			lowPassFilterMap(Mrefs[iref], lowpass, angpix);
		}
	}

	if (fn_mask != "")
	{
		for(int iref=0; iref < Mrefs.size();iref++)
		{
			Image<DOUBLE> msk;
			msk.read(fn_mask);
			msk().setXmippOrigin();
			if (!msk().sameShape(Mrefs[iref]))
				REPORT_ERROR("project ERROR: mask and map have different sizes!");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mrefs[iref])
				DIRECT_MULTIDIM_ELEM(Mrefs[iref], n) *= DIRECT_MULTIDIM_ELEM(msk(), n);
		}
	}

	if (do_invert_refs)
	{
		for(int iref=0; iref < Mrefs.size();iref++)
		{
			Mrefs[iref] *= -1;
		}
	}
	size_refs.clear();
	for(int iref=0; iref < Mrefs.size();iref++)
	{
		size_refs.push_back(XSIZE(Mrefs[iref]));
	}
	PPref.clear();
	MultidimArray<DOUBLE> dummy;
	for(int iref=0; iref < Mrefs.size();iref++)
	{
		Projector projector(size_refs[iref]);
		projector.computeFourierTransformMap(Mrefs[iref],dummy,size_refs[iref]);
		PPref.push_back(projector);
	}

	// background radius of new rewindowed particles
	bg_radius = (int) (new_particle_diameter / 2. / angpix); 
}


void SubtractingMicrograph::putAllProjectsIntoEmptyMicrograph(Image<DOUBLE> &emptyMic,MetaDataTable &MD,int Maxdim,int Xdim, int Ydim,DOUBLE binning,bool do_apply_ctf_to_total)
{
	FourierTransformer transformer;
	MultidimArray<DOUBLE> Mref_rot, Fctf;
	MultidimArray<Complex> Fref;
	DOUBLE psi,tilt,rot,xcoord,ycoord;
	int iref,x_empty,y_empty;
	Matrix1D<DOUBLE> offset(2);
	Matrix1D<DOUBLE> p_ptcls,p_ref;
	Matrix2D<DOUBLE> A;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		// For Euler angles 
		MD.getValue(EMDL_ORIENT_ROT,rot);
		MD.getValue(EMDL_ORIENT_TILT,tilt);
		MD.getValue(EMDL_ORIENT_PSI,psi);
		Euler_angles2matrix(rot,tilt,psi,A);
		// find offset of particle
		MD.getValue(EMDL_ORIENT_ORIGIN_X,XX(offset));
		MD.getValue(EMDL_ORIENT_ORIGIN_Y,YY(offset));
		// Get coordinates
		MD.getValue(EMDL_IMAGE_COORD_X,xcoord);
		MD.getValue(EMDL_IMAGE_COORD_Y,ycoord);
		// For class number
		if (!MD.getValue(EMDL_PARTICLE_CLASS, iref) && PPref.size() > 1)
	                REPORT_ERROR("subtract_micrograph: cannot find class number in input STAR file");
		iref -= 1 ;
		Fref.resize(size_refs[iref],size_refs[iref] / 2 + 1);
		PPref[iref].get2DFourierTransform(Fref,A,IS_NOT_INV);
		// apply offset of particle to reference
		if (ABS(XX(offset)) > 0. || ABS(YY(offset)) > 0. )
		{
			shiftImageInFourierTransform(Fref,Fref,tab_sin,tab_cos,(DOUBLE)size_refs[iref],-XX(offset),-YY(offset));
		}
		// apply ctf to projections if not apply it to micrograph
		if (do_ctf && (! do_apply_ctf_to_total) )
		{
			CTF ctf;
			ctf.read(MD, MD);
			Fctf.resize(Fref);
			ctf.getFftwImage(Fctf, size_refs[iref], size_refs[iref], angpix, false, false, intact_ctf_first_peak, true);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
			{
				DIRECT_MULTIDIM_ELEM(Fref, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
			}
		}
		// invert to real space
		Mref_rot.resize(size_refs[iref],size_refs[iref]);
		transformer.inverseFourierTransform(Fref,Mref_rot);
		CenterFFT(Mref_rot, true);
		//Mref_rot.setXmippOrigin();
		// Get the integer coordinates in emptyMic from that in raw micrograph
		x_empty = (int)(xcoord / binning ) + CEIL((DOUBLE)( Maxdim - Xdim) / 2.)  ;
		y_empty = (int)(ycoord / binning ) + CEIL((DOUBLE)( Maxdim - Ydim) / 2.)  ;
		// put projection into empty micrograph
		// i=0...Ysize,j=0...Xsize
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mref_rot)
		{
			if ((j + x_empty - size_refs[iref] / 2 ) >= 0 && ( j + x_empty - size_refs[iref] / 2 ) < Maxdim && (i + y_empty - size_refs[iref] / 2) >= 0 && ( i + y_empty - size_refs[iref] / 2) < Maxdim )
			{
				DIRECT_A2D_ELEM(emptyMic(),i + y_empty - size_refs[iref] / 2,j + x_empty - size_refs[iref] / 2) += DIRECT_A2D_ELEM(Mref_rot,i,j);
			}
		}
	}
}

void SubtractingMicrograph::updateCoord2Target(MetaDataTable &MD,DOUBLE binning)
{
	DOUBLE psi,tilt,rot,xcoord,ycoord,newxcoord,newycoord;
	Matrix1D<DOUBLE> offset(2);
	Matrix1D<DOUBLE> p_ptcls,p_ref;
	Matrix2D<DOUBLE> A;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		// For Euler angles 
		MD.getValue(EMDL_ORIENT_ROT,rot);
		MD.getValue(EMDL_ORIENT_TILT,tilt);
		MD.getValue(EMDL_ORIENT_PSI,psi);
		Euler_angles2matrix(rot,tilt,psi,A);
		// find offset of particle
		MD.getValue(EMDL_ORIENT_ORIGIN_X,XX(offset));
		MD.getValue(EMDL_ORIENT_ORIGIN_Y,YY(offset));
		// Get coordinates
		MD.getValue(EMDL_IMAGE_COORD_X,xcoord);
		MD.getValue(EMDL_IMAGE_COORD_Y,ycoord);
		// Update the coordinates and origin shifts of target

		p_ref = vectorR3(target_x,target_y,target_z) ;
		p_ptcls = A * p_ref ;
		// new coord in rebinningd micrographs
		newxcoord = floor(xcoord  - XX(offset) * binning + XX(p_ptcls) * binning );
		newycoord = floor(ycoord  - YY(offset) * binning + YY(p_ptcls) * binning );
		MD.setValue(EMDL_ORIENT_ORIGIN_X,0.);
		MD.setValue(EMDL_ORIENT_ORIGIN_Y,0.);
		MD.setValue(EMDL_IMAGE_COORD_X,newxcoord);
		MD.setValue(EMDL_IMAGE_COORD_Y,newycoord);
	}
}



void SubtractingMicrograph::rewindowRecenteredParticles(Image<DOUBLE> &Imic,const int icmic)
{
	int ipos = 0;
	long int total_ipos;
	DOUBLE all_avg,all_stddev,all_minval,all_maxval;
	Image<DOUBLE> Ipart,neighbor_Ipart,emptyNeighborParticle;

	total_ipos=MD_coord_mics[icmic].numberOfObjects();
	FileName fn_output_img;
	fn_output_img="Particles/" + fn_coord_mics[icmic].withoutExtension() + "_" + fn_rewindow + ".mrcs";

	// Copy Current MetaDataTable
	MetaDataTable MD_copy;
	MD_copy=MD_coord_mics[icmic];
	//long int dummy;
	//dummy = MD_coord_mics[icmic].firstObject();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_coord_mics[icmic])
	{
		DOUBLE dxpos, dypos;
		long int xpos, ypos;
		long int x0, xF, y0, yF;
		MD_coord_mics[icmic].getValue(EMDL_IMAGE_COORD_X, dxpos);
		MD_coord_mics[icmic].getValue(EMDL_IMAGE_COORD_Y, dypos);
		// xpos, ypos: integer with reference binning factor
		xpos = (long int) (dxpos / binning);
		ypos = (long int) (dypos / binning);
		x0 = xpos  +  FIRST_XMIPP_INDEX(extract_size);
		xF = xpos  +  LAST_XMIPP_INDEX(extract_size);
		y0 = ypos  +  FIRST_XMIPP_INDEX(extract_size);
		yF = ypos  +  LAST_XMIPP_INDEX(extract_size);
		// Do not discard particles that are completely outside the micrograph and print a warning
		// To simplify the matching between the number of the original and rewindowed particles
		if (yF < 0 || y0 >= YSIZE(Imic()) || xF < 0 || x0 >= XSIZE(Imic()) )
		{
			std::cout << "Warning: particle" + integerToString(ipos+1) + "@" + fn_output_img +  " lies completely outside micrograph. Fill it with the micrograph center" << std::endl ;
			xpos = XSIZE(Imic()) / 2;
			ypos = YSIZE(Imic()) / 2;
			MD_coord_mics[icmic].setValue(EMDL_IMAGE_COORD_X, xpos * binning );
			MD_coord_mics[icmic].setValue(EMDL_IMAGE_COORD_Y, ypos * binning );
			x0 = xpos  +  FIRST_XMIPP_INDEX(extract_size);
			xF = xpos  +  LAST_XMIPP_INDEX(extract_size);
			y0 = ypos  +  FIRST_XMIPP_INDEX(extract_size);
			yF = ypos  +  LAST_XMIPP_INDEX(extract_size);
		}

		if ( do_subtract_neighbor)
		{
			x0 = xpos +  FIRST_XMIPP_INDEX(neighbor_size);
			xF = xpos +  LAST_XMIPP_INDEX(neighbor_size);
			y0 = ypos +  FIRST_XMIPP_INDEX(neighbor_size);
			yF = ypos +  LAST_XMIPP_INDEX(neighbor_size) ;
			// extract one padding particle
			Imic().window(neighbor_Ipart(), y0, x0, yF, xF);

			MetaDataTable MD_neighborParticles;
			// Find neighbor particles
			int neighbor_ipos = 0;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_copy)
			{
				if (neighbor_ipos != ipos)
				{
					DOUBLE Ndxpos, Ndypos;
					long int Nxpos, Nypos,dx,dy,ndist;
					MD_copy.getValue(EMDL_IMAGE_COORD_X, Ndxpos);
					MD_copy.getValue(EMDL_IMAGE_COORD_Y, Ndypos);
					// Nxpos, Nypos: integer with reference binning factor
					Nxpos = (long int)(Ndxpos / binning);
					Nypos = (long int)(Ndypos / binning);

					int iref;
					MD_copy.getValue(EMDL_PARTICLE_CLASS, iref);
					iref -= 1;
					ndist = (neighbor_size + size_refs[iref] ) / 2 ;
					dx=Nxpos  - xpos;
					dy=Nypos  - ypos;
					if ( dx < ndist && dx > -ndist && dy < ndist && dy > -ndist )
					{
						MD_neighborParticles.addObject(MD_copy.getObject());
						MD_neighborParticles.setValue(EMDL_IMAGE_COORD_X, Ndxpos - dxpos + (neighbor_size / 2 * binning));
						MD_neighborParticles.setValue(EMDL_IMAGE_COORD_Y, Ndypos - dypos + (neighbor_size / 2 * binning));
					}
				}
				neighbor_ipos++;
			}
	
			emptyNeighborParticle().resize(neighbor_size,neighbor_size);
			emptyNeighborParticle().initZeros();
			//dummy = MD_neighborParticles.firstObject();
			putAllProjectsIntoEmptyMicrograph(emptyNeighborParticle,MD_neighborParticles,neighbor_size,neighbor_size,neighbor_size,binning);
			// Apply ctf to micrograph, assuming no defocus varying 
			if (do_ctf)
			{
				FourierTransformer transformer;
				MultidimArray<Complex> FemptyNeighborParticle;
				FemptyNeighborParticle.resize(neighbor_size, neighbor_size / 2 + 1 );
				transformer.FourierTransform(emptyNeighborParticle(),FemptyNeighborParticle);
				MultidimArray<DOUBLE> Fctf;
				Fctf.resize(FemptyNeighborParticle);
				ctf_mics[icmic].getFftwImage(Fctf, neighbor_size, neighbor_size, angpix, false, false, intact_ctf_first_peak, true);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FemptyNeighborParticle)
				{
					DIRECT_MULTIDIM_ELEM(FemptyNeighborParticle, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
				transformer.inverseFourierTransform(FemptyNeighborParticle,emptyNeighborParticle());
			}
			//emptyNeighborParticle().setXmippOrigin();
			//windowed images have Xmipp origin.
			STARTINGX(neighbor_Ipart()) = 0;
			STARTINGY(neighbor_Ipart()) = 0;
		
			neighbor_Ipart() -= emptyNeighborParticle();
			emptyNeighborParticle.clear();

			// extract one particle in Ipart
			x0 = neighbor_size / 2 + FIRST_XMIPP_INDEX(extract_size);
			xF = neighbor_size / 2 + LAST_XMIPP_INDEX(extract_size);
			y0 = neighbor_size / 2 + FIRST_XMIPP_INDEX(extract_size);
			yF = neighbor_size / 2 + LAST_XMIPP_INDEX(extract_size);
			neighbor_Ipart().window(Ipart(), y0, x0, yF, xF);

			x0 = xpos  +  FIRST_XMIPP_INDEX(extract_size);
			xF = xpos  +  LAST_XMIPP_INDEX(extract_size);
			y0 = ypos  +  FIRST_XMIPP_INDEX(extract_size);
			yF = ypos  +  LAST_XMIPP_INDEX(extract_size);
		}
		else
		{
			// extract one particle in Ipart
			Imic().window(Ipart(), y0, x0, yF, xF);
		}
		// Check boundaries: fill pixels outside the boundary with the nearest ones inside
		// This will create lines at the edges, rather than zeros
		Ipart().setXmippOrigin();

		// X-boundaries
		if (x0 < 0 || xF >= XSIZE(Imic()) )
		{
			FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart())
			{
				if (j + xpos < 0)
					A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, i, -xpos);
				if (j + xpos >= XSIZE(Imic()))
					A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, i, XSIZE(Imic()) - xpos - 1);
			}
		}

		// Y-boundaries
		if (y0 < 0 || yF >= YSIZE(Imic()))
		{
			FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart())
			{
				if (i + ypos < 0)
					A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, -ypos, j);
				if (i + ypos >= YSIZE(Imic()))
					A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, YSIZE(Imic()) - ypos - 1, j);
			}
		}

		//
		// performPerImageOperations will also append the particle to the output stack in fn_stack
		Ipart().setXmippOrigin();
	
		if (do_normalise) normalise(Ipart, bg_radius, white_dust_stddev, black_dust_stddev,do_ramp);

		//if (do_invert_contrast) invert_contrast(Ipart);

		// Calculate mean, stddev, min and max
		DOUBLE avg, stddev, minval, maxval;
		Ipart().computeStats(avg, stddev, minval, maxval);


		// Keep track of overall statistics
		all_minval = XMIPP_MIN(minval, all_minval);
		all_maxval = XMIPP_MAX(maxval, all_maxval);
		all_avg	+= avg;
		all_stddev += stddev*stddev;

		// Last particle: reset the min, max, avg and stddev values in the main header
		if (ipos == (total_ipos - 1))
		{
			all_avg /= total_ipos;
			all_stddev = sqrt(all_stddev/total_ipos);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, all_minval);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, all_maxval);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, all_avg);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, all_stddev);
		}

		// Write this particle to the stack on disc
		// First particle: write stack in overwrite mode, from then on just append to it
		if (ipos == 0)
			Ipart.write(fn_output_img, -1, false, WRITE_OVERWRITE);
		else
			Ipart.write(fn_output_img, -1, false, WRITE_APPEND);
		// Update image name of rewindowed particles
		FileName fn_rewindowimg;
		fn_rewindowimg.compose(ipos + 1,fn_output_img);
		MD_coord_mics[icmic].setValue(EMDL_IMAGE_NAME,fn_rewindowimg);
		ipos++;
	}

	FileName fn_star;
	fn_star =  "Particles/" + fn_coord_mics[icmic].withoutExtension() + "_" + fn_rewindow + ".star";
	MD_coord_mics[icmic].write(fn_star);
}

void SubtractingMicrograph::subtractProjectionsFromParticles(const int imic)
{
	Image<DOUBLE> Ipart,emptyPart,Ipartsub;
	MetaDataTable MD_newmics,MD_newsub_mics,MD_newcoord_mics;
	// image name
	FileName fn_img,fn_sub_img,fn_coord_img,fn_stack;
	// output stack name
	FileName fn_img_out, fn_sub_img_out, fn_coord_img_out;
	fn_img_out="Particles/" + fn_mics[imic].withoutExtension() + "_raw_" + fn_rewindow + ".mrcs";
	fn_sub_img_out="Particles/" + fn_mics[imic].withoutExtension() + "_sub_" + fn_rewindow + ".mrcs";
	fn_coord_img_out="Particles/" + fn_mics[imic].withoutExtension() + "_coord_" + fn_rewindow + ".mrcs";
	long int num_in,num_sub,num_coord;
	int ipos;
	long int total_ipos;
	DOUBLE all_avg,all_stddev,all_minval,all_maxval;
	DOUBLE Ndxpos, Ndypos, dxpos, dypos;

	int coord_imic =  -1;
	for (int icmic = 0; icmic < fn_coord_mics.size(); icmic++)
	{
		if (fn_coord_mics[icmic] == fn_mics[imic])
		{
			coord_imic=icmic;
			break;
		}
	}
	int sub_imic =  -1;
	for (int icmic = 0; icmic < fn_sub_mics.size(); icmic++)
	{
		if (fn_sub_mics[icmic] == fn_mics[imic])
		{
			sub_imic=icmic;
			break;
		}
	}
	if ( coord_imic >= 0 && sub_imic >= 0 )
	{
		//Filling three new metadata_tables
		MD_newmics.clear();
		MD_newsub_mics.clear();
		MD_newcoord_mics.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_mics[imic])
		{
			MD_mics[imic].getValue(EMDL_IMAGE_NAME, fn_img);
			fn_img.decompose(num_in,fn_stack);
	
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_sub_mics[sub_imic])
			{
				MD_sub_mics[sub_imic].getValue(EMDL_IMAGE_NAME, fn_sub_img);
				fn_sub_img.decompose(num_sub,fn_stack);
				
				if (num_sub == num_in)
				{
					FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_coord_mics[coord_imic])
					{
						MD_coord_mics[coord_imic].getValue(EMDL_IMAGE_NAME, fn_coord_img);
						fn_coord_img.decompose(num_coord,fn_stack);
				
						if (num_coord == num_in)
						{
							MD_newmics.addObject(MD_mics[imic].getObject());
							MD_newsub_mics.addObject(MD_sub_mics[sub_imic].getObject());
							MD_newcoord_mics.addObject(MD_coord_mics[coord_imic].getObject());
							break;
						}
					}
					break;
				}
			}
		}
		total_ipos = MD_newmics.numberOfObjects();
		ipos = 0;
		updateCoord2Target(MD_newcoord_mics,binning);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_newmics)
		{
			MD_newmics.getValue(EMDL_IMAGE_NAME, fn_img);
			Ipart.read(fn_img);
			if (do_output_raw_particles)
			{
				fn_img.compose(ipos + 1,fn_img_out);
				MD_newmics.setValue(EMDL_IMAGE_NAME, fn_img);
				if (ipos == 0)
					Ipart.write(fn_img_out, -1, false, WRITE_OVERWRITE);
				else
					Ipart.write(fn_img_out, -1, false, WRITE_APPEND);
			}

			STARTINGX(Ipart()) =  0;
			STARTINGY(Ipart()) =  0;
			MD_newmics.getValue(EMDL_IMAGE_COORD_X, dxpos);
			MD_newmics.getValue(EMDL_IMAGE_COORD_Y, dypos);
			// subtract projection from Raw image
			if ( do_subtract )
			{
				MetaDataTable MD_one;
				MD_one.addObject(MD_newsub_mics.getObject(ipos));
				
				if (do_output_sub_particles)
				{
					MD_newsub_mics.getValue(EMDL_IMAGE_NAME, fn_sub_img,ipos);
					Ipartsub.read(fn_sub_img);
					fn_sub_img.compose(ipos + 1,fn_sub_img_out);
					MD_newsub_mics.setValue(EMDL_IMAGE_NAME, fn_sub_img,ipos);
					if (ipos == 0)
						Ipartsub.write(fn_sub_img_out, -1, false, WRITE_OVERWRITE);
					else
						Ipartsub.write(fn_sub_img_out, -1, false, WRITE_APPEND);
				}
				long int Xdim,Ydim,Zdim,Ndim;
				Ipart().getDimensions(Xdim,Ydim,Zdim,Ndim);

				MD_one.getValue(EMDL_IMAGE_COORD_X, Ndxpos,0);
				MD_one.getValue(EMDL_IMAGE_COORD_Y, Ndypos,0);
				
				MD_one.setValue(EMDL_IMAGE_COORD_X, floor((Ndxpos - dxpos) / binning ) + (Xdim / 2. ),0);
				MD_one.setValue(EMDL_IMAGE_COORD_Y, floor((Ndypos - dypos) / binning ) + (Ydim / 2. ),0);
	
				emptyPart().resize(Xdim,Ydim);
				emptyPart().initZeros();
				putAllProjectsIntoEmptyMicrograph(emptyPart,MD_one,Xdim,Xdim,Ydim,binning);
				// Apply ctf to micrograph, assuming no defocus varying 
				if (do_ctf)
				{
					CTF ctf_one;
					ctf_one.read(MD_one,MD_one,0);
					FourierTransformer transformer;
					MultidimArray<Complex> FemptyPart;
					FemptyPart.resize(Xdim, Xdim / 2 + 1 );
					transformer.FourierTransform(emptyPart(),FemptyPart);
					MultidimArray<DOUBLE> Fctf;
					Fctf.resize(FemptyPart);
					ctf_one.getFftwImage(Fctf, Xdim, Xdim, angpix, false, false, intact_ctf_first_peak, true);
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FemptyPart)
					{
						DIRECT_MULTIDIM_ELEM(FemptyPart, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
					transformer.inverseFourierTransform(FemptyPart,emptyPart());
				}
				Ipart() -= emptyPart();
			}
			Ipart().setXmippOrigin();
					
			MD_newcoord_mics.getValue(EMDL_IMAGE_COORD_X, Ndxpos,ipos);
			MD_newcoord_mics.getValue(EMDL_IMAGE_COORD_Y, Ndypos,ipos);
	
			Matrix1D<DOUBLE> refsuboffset(2);
			XX(refsuboffset) = floor((dxpos - Ndxpos) / binning );
			YY(refsuboffset) = floor((dypos - Ndypos) / binning );

			MD_newcoord_mics.setValue(EMDL_IMAGE_COORD_X, dxpos - (XX(refsuboffset) * binning),ipos);
			MD_newcoord_mics.setValue(EMDL_IMAGE_COORD_Y, dypos - (YY(refsuboffset) * binning),ipos);
			selfTranslate(Ipart(),refsuboffset,true);
			rewindow(Ipart,extract_size);
		
			if (do_normalise) normalise(Ipart, bg_radius, white_dust_stddev, black_dust_stddev,do_ramp);
	
			// Calculate mean, stddev, min and max
			DOUBLE avg, stddev, minval, maxval;
			Ipart().computeStats(avg, stddev, minval, maxval);
			// Keep track of overall statistics
			all_minval = XMIPP_MIN(minval, all_minval);
			all_maxval = XMIPP_MAX(maxval, all_maxval);
			all_avg	+= avg;
			all_stddev += stddev*stddev;
	
			// Last particle: reset the min, max, avg and stddev values in the main header
			if (ipos == (total_ipos - 1))
			{
				all_avg /= total_ipos;
				all_stddev = sqrt(all_stddev/total_ipos);
				Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, all_minval);
				Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, all_maxval);
				Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, all_avg);
				Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, all_stddev);
			}
	
			// Write this particle to the stack on disc
			// First particle: write stack in overwrite mode, from then on just append to it
			fn_coord_img.compose(ipos + 1,fn_coord_img_out);
			MD_newcoord_mics.setValue(EMDL_IMAGE_NAME,fn_coord_img,ipos);
			if (ipos == 0)
				Ipart.write(fn_coord_img_out, -1, false, WRITE_OVERWRITE);
			else
				Ipart.write(fn_coord_img_out, -1, false, WRITE_APPEND);
			ipos++;
		}
	
		FileName fn_star;
		if (do_output_raw_particles)
		{
			fn_star =  "Particles/" + fn_mics[imic].withoutExtension() + "_raw_" + fn_rewindow + ".star";
			MD_newmics.write(fn_star);
		}
		if (do_subtract && do_output_sub_particles)
		{
			fn_star =  "Particles/" + fn_sub_mics[sub_imic].withoutExtension() + "_sub_" + fn_rewindow + ".star";
			MD_newsub_mics.write(fn_star);
		}
		fn_star =  "Particles/" + fn_coord_mics[coord_imic].withoutExtension() + "_coord_" + fn_rewindow + ".star";
		MD_newcoord_mics.write(fn_star);
	}
}

void SubtractingMicrograph::subtractProjectionsFromMicrograph(const int imic)
{
	Image<DOUBLE> Imic,emptyMic,centor_emptyMic;
	Imic.read(fn_mics[imic]);
	long int Xdim,Ydim,Zdim;
	long int Ndim;
	Imic().getDimensions(Xdim,Ydim,Zdim,Ndim);

	if (ABS(binning - 1.0) > 0.001)
	{
		long int intbinningd_x,intbinningd_y,int_x,int_y;
		intbinningd_x = (long int) floor((DOUBLE)Xdim / binning) ;
		intbinningd_y = (long int) floor((DOUBLE)Ydim / binning) ;
		int_x = (long int) ROUND((DOUBLE)intbinningd_x * binning);
		int_y = (long int) ROUND((DOUBLE)intbinningd_y * binning);
		if ((Xdim - int_x) > 0 || (Ydim - int_y) > 0 )
		{
			MultidimArray<DOUBLE> intTrunc;
			intTrunc.resize(int_y,int_x);
			Imic().window(intTrunc,0,0,int_y,int_x);
			selfScaleToSize(intTrunc, intbinningd_x, intbinningd_y);
			Imic().resize(intbinningd_y,intbinningd_x);
			Imic()=intTrunc;
			Imic().getDimensions(Xdim,Ydim,Zdim,Ndim);
		}
		else
		{
			selfScaleToSize(Imic(), (int)((DOUBLE)Xdim / binning), (int)((DOUBLE)Ydim / binning));
			Imic().getDimensions(Xdim,Ydim,Zdim,Ndim);
		}
	}

	// normalise micrograph
	normalise(Imic, XMIPP_MIN(Xdim,Ydim) / 2, white_dust_stddev, black_dust_stddev,do_ramp);
	if (do_invert_micrographs) invert_contrast(Imic);
	
	if ( do_subtract)
	{
		int Maxdim = XMIPP_MAX(Xdim,Ydim);
		// for coherent scaling between micrographs and particles
		if (Maxdim % 2 != 0) Maxdim += 1;
		emptyMic().resize(Maxdim,Maxdim);
		emptyMic().initZeros();
		putAllProjectsIntoEmptyMicrograph(emptyMic,MD_mics[imic],Maxdim,Xdim,Ydim,binning,do_apply_ctf_to_total);
		// Apply ctf to micrograph, assuming no defocus varying 
		if (do_ctf && do_apply_ctf_to_total)
		{
			FourierTransformer transformer;
			MultidimArray<Complex> FemptyMic;
			FemptyMic.resize(Maxdim, Maxdim / 2 + 1 );
			transformer.FourierTransform(emptyMic(),FemptyMic);
			MultidimArray<DOUBLE> Fctf;
			Fctf.resize(FemptyMic);
			ctf_mics[imic].getFftwImage(Fctf, Maxdim, Maxdim, angpix, false, false, intact_ctf_first_peak, true);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FemptyMic)
			{
				DIRECT_MULTIDIM_ELEM(FemptyMic, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
			}
			transformer.inverseFourierTransform(FemptyMic,emptyMic());
		}
		long int x0 = Maxdim / 2 + FIRST_XMIPP_INDEX(Xdim);
		long int xF = Maxdim / 2 + LAST_XMIPP_INDEX(Xdim);
		long int y0 = Maxdim / 2 + FIRST_XMIPP_INDEX(Ydim);
		long int yF = Maxdim / 2 + LAST_XMIPP_INDEX(Ydim);
		// extract centor of emptyMic with size of Imic
		emptyMic().window(centor_emptyMic(), y0, x0, yF, xF);
	
		STARTINGX(centor_emptyMic()) = 0;
		STARTINGY(centor_emptyMic()) = 0;
		if (do_output_projection_micrographs)
		{
			FileName fn_output_projmic;
			fn_output_projmic=fn_mics[imic].insertBeforeExtension("_" + fn_projectionmic);
			centor_emptyMic.write(fn_output_projmic);
		}
		Imic() -= centor_emptyMic();
		emptyMic.clear();
		centor_emptyMic.clear();
		if (do_output_refsub_micrographs)
		{
			FileName fn_output_mic;
			fn_output_mic=fn_mics[imic].insertBeforeExtension("_" + fn_refsubmic);
			Imic.write(fn_output_mic);
		}
	}

	int coord_imic =  -1;
	for (int icmic = 0; icmic < fn_coord_mics.size(); icmic++)
	{
		if (fn_coord_mics[icmic] == fn_mics[imic])
		{
			coord_imic=icmic;
			break;
		}
	}
	if (coord_imic >= 0)
	{
		if (do_reextract || do_subtract_neighbor)
		{
			updateCoord2Target(MD_coord_mics[coord_imic],binning);
			rewindowRecenteredParticles(Imic,coord_imic);
		}
	}
}



void SubtractingMicrograph::joinStarFiles( std::vector<FileName> &fn_mics, FileName fn_suffix, FileName fn_output)
{
	MetaDataTable MDout, MDone;
	FileName fn_star;
	for (int imic=0; imic < fn_mics.size();imic++)
	{
		fn_star = "Particles/" + fn_mics[imic].withoutExtension() + fn_suffix;
		if (exists(fn_star))
		{
			MDone.read(fn_star);
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDone)
			{
				MDout.addObject(MDone.getObject());
			}
		}
	}
	MDout.write(fn_output);
	std::cout << " Written out STAR file with all particles in " << fn_output << std::endl;
}

void SubtractingMicrograph::renumberInputStarFileAndOutputRawImages(const int imic)
{
	FileName fn_img,fn_stack,fn_img_renumber;
	long int dummy;
	Image<DOUBLE> Iimg;
	int n_img = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_mics[imic])
	{
		MD_mics[imic].getValue(EMDL_IMAGE_NAME, fn_img);
		Iimg.read(fn_img);
		fn_img.decompose(dummy,fn_stack);
		fn_img_renumber.compose(n_img + 1,fn_stack.insertBeforeExtension("_" + fn_renumber));
		if (n_img == 0)
			Iimg.write(fn_img_renumber,-1,false,WRITE_OVERWRITE);
		else
			Iimg.write(fn_img_renumber,-1,false,WRITE_APPEND);
		MD_mics[imic].setValue(EMDL_IMAGE_NAME, fn_img_renumber);
		n_img++;
	}
	FileName fn_star;
	fn_star = "Particles/" + fn_mics[imic].withoutExtension() + "_" + fn_renumber + ".star";
	MD_mics[imic].write(fn_star);
}


void SubtractingMicrograph::run()
{
	int nr_mics=fn_mics.size();
	
	int barstep;
	if ( do_renumberRawImages )
	{
		if (verb > 0 && do_renumberRawImages )
		{
			std::cout << "Renumbering all particles in the input star file..." << std::endl;
			init_progress_bar(nr_mics);
			barstep = XMIPP_MAX(1, nr_mics / 60);
		}
		for (int imic=0; imic < fn_mics.size();imic++)
		{
			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
			renumberInputStarFileAndOutputRawImages(imic);
		}
		joinStarFiles(fn_mics, "_" + fn_renumber + ".star",fn_in.insertBeforeExtension("_" + fn_renumber) );
	}
	if (do_rewindowMicrographs)
	{
		if (verb > 0)
		{
			std::cout << "Processing all micrographs in the input star file..." << std::endl;
			init_progress_bar(nr_mics);
			barstep = XMIPP_MAX(1, nr_mics / 60);
		}
		for (int imic=0; imic < fn_mics.size();imic++)
		{
			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
			subtractProjectionsFromMicrograph(imic);
		}
		if (do_reextract || do_subtract_neighbor)
		{
			joinStarFiles(fn_mics,"_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_" + fn_rewindow) );
		}
	}

	if (do_rewindowParticles)
	{
		if (verb > 0)
		{
			std::cout << "Processing all particles in the input star file..." << std::endl;
			init_progress_bar(nr_mics);
			barstep = XMIPP_MAX(1, nr_mics / 60);
		}
		for (int imic=0; imic < fn_mics.size();imic++)
		{
			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
			subtractProjectionsFromParticles(imic);
		}
		if (do_output_raw_particles)
		{
			joinStarFiles(fn_mics, "_raw_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_raw_" + fn_rewindow) );
		}
		if (do_output_sub_particles)
		{
			joinStarFiles(fn_mics, "_sub_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_sub_" + fn_rewindow) );
		}
		joinStarFiles(fn_mics, "_coord_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_coord_" + fn_rewindow) );
	}

}

