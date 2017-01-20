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


#ifndef SUBTRACTINGMICROGRAPH_H_
#define SUBTRACTINGMICROGRAPH_H_
#include  <glob.h>
#include  <vector>
#include  <string>
#include  <stdlib.h>
#include  <stdio.h>
#include "src/ctf.h"
#include "src/image.h"
#include "src/projector.h"
#include "src/complex.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/ctffind_runner.h"
#include "src/parallel.h"
#include <src/fftw.h>
#include <src/time.h>

class SubtractingMicrograph
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input file name
	FileName fn_in;

	// Input ref file name
	FileName fn_ref;

	// Input mask file name
	FileName fn_mask;

	// Input subtract file name
	FileName fn_sub;

	// Input coordinate file name
	FileName fn_coord;

	// MetadataTable storing input star file
	MetaDataTable MDin,MD_coord,MD_sub;

	// Output rootname
	FileName fn_renumber, fn_rewindow,fn_refsubmic,fn_projectionmic;

	// Do renumber raw images?
	bool do_renumberRawImages;

	// only rewindow raw images but not subtract projections?
	bool do_onlyRewindow;

	// Pixel size of micrographs
	DOUBLE angpix;

	// low pass filter reference map
	DOUBLE lowpass;

	//binning factor of references relative to micrographs
	DOUBLE binning;

	// apply ctf to projections of references?
	bool do_ctf;

	// apply ctf to sum micrograph containing total projections?
	bool do_apply_ctf_to_total;
	
	// intact ctf first peak
	bool intact_ctf_first_peak;

	// References have inverted contrast?
	bool do_invert_refs;
	// Micrographs have inverted contrast?
	bool do_invert_micrographs;

	// Output projection-subtracted micrographs?
	bool do_output_refsub_micrographs;

	// Output projection micrographs?
	bool do_output_projection_micrographs;

	// Do not re-output raw particles
	bool do_output_raw_particles;

	// Do not re-output subtraction particles
	bool do_output_sub_particles;

	// Do rewindow micrographs?
	bool do_rewindowMicrographs;

	// Do rewindow particles?
	bool do_rewindowParticles;

	// Subtract projection before re-extraction?
	bool do_subtract;

	// Perform particle re-extraction?
	bool do_reextract;

	// Box size to extract the particles in
	int extract_size;

	// Radius of a circle in the extracted images outside of which one calculates background mean and stddev
	int bg_radius;

	// Particle diameter of target part to be rewindowed
	DOUBLE new_particle_diameter;

	// coordinate X,Y,Z of target in original reference map
	// The origin is at the centor of the map.
	DOUBLE target_x,target_y,target_z;

	// Bias in picked coordinates in X and in Y direction (in pixels)
	//DOUBLE extract_bias_x, extract_bias_y;

	////////////////////////////////////// Post-extraction image modifications
	// Perform re-scaling of extracted images
	//bool do_rescale;

	// Perform re-windowing of extracted images
	//bool do_rewindow;
	//int window;

	// Perform normalization of the extract images
	bool do_normalise;

	//Subtract ramp instead of a level background in normalization
	bool do_ramp;

	// Standard deviations to remove black and white dust
	DOUBLE white_dust_stddev, black_dust_stddev;

	//////////////////////////////////// Output STAR file
	bool do_join_starfile;

	// Tabulated sin and cosine functions for shifts in Fourier space
	TabSine tab_sin;
	TabCosine tab_cos;

	// store all of micrographs name
	std::vector<FileName> fn_mics,fn_coord_mics,fn_sub_mics,fn_allmics,fn_mymics;

	// CTF models of all micrographs
	std::vector<CTF> ctf_mics,ctf_coord_mics,ctf_sub_mics; 

	// Split MD into fractions
	std::vector<MetaDataTable> MD_mics, MD_coord_mics,MD_sub_mics;

	// all references
	std::vector<MultidimArray<DOUBLE> > Mrefs;

	// Fourier transform of references
	std::vector<Projector> PPref;

	// The size of references
	std::vector<int> size_refs;

	//Whether subtract neighbor particles
	bool do_subtract_neighbor;

	// The size of neighbor
	int neighbor_size;

public:
	// Read command line arguments
	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise(int rank = 0, int nodesize = 0 );

	// General Running
	void run();

	// Put projections of references into a temporary image container
	// The Micrograph name, coordinates and origin offset of particles 
	// are updated in this step.
	void putAllProjectsIntoEmptyMicrograph(Image<DOUBLE> &emptyMic,MetaDataTable &MD,int Maxdim,int Xdim, int Ydim,DOUBLE binning,bool do_apply_ctf_to_total=true);

	// update coordinates to target part
	void updateCoord2Target(MetaDataTable &MD,DOUBLE binning);

	// Rewindow particles from reference projection-subtracted micrographs
	// The image name of particles are updated in this step.
	void rewindowRecenteredParticles(Image<DOUBLE> &Imic,const int icmic);

	// Processing each micrograph
	void subtractProjectionsFromMicrograph(const int imic);

	// Processing each micrograph
	void subtractProjectionsFromParticles(const int imic);

	// Join all star files
	void joinStarFiles(std::vector<FileName> &fn_mics, FileName fn_suffix, FileName fn_output);

	// Renumber original particles and output them
	void renumberInputStarFileAndOutputRawImages(const int imic);
	
	// Split star file
	void splitMetaDataAccordingtoMicrograph(MetaDataTable &MD, std::vector<MetaDataTable> &MD_mics, std::vector<FileName> &fn_mics, std::vector<CTF> &ctf_mics, std::vector<FileName> &fn_mymics);
	// Get all micrographs 
	void getMicrographsFromMetaData(MetaDataTable &MD, std::vector<FileName> &fn_allmics );
};

#endif /* SUBTRACTINGMICROGRAPH_H_ */



