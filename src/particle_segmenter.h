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

#ifndef PARTICLE_SEGMENTER_H_
#define PARTICLE_SEGMENTER_H_

#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/projector.h"
#include "src/ctf.h"
#include <src/fftw.h>
#include <src/time.h>
#include <src/funcs.h>


#define WIDTH_FMASK_EDGEB 2

#define REWINDOWOFFSETX  0
#define REWINDOWOFFSETY  1
#define NR_REWINDOWOFFSETS 2


class ParticleSegmenter
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input & Output rootname
	FileName fn_in, fn_ref, fn_out;

	// Input metadata
	MetaDataTable MDin;

	// Pixel size (for low-pass filter and particle diameter)
	DOUBLE angpix;

	// Particle diameter (in Angstroms)
	DOUBLE particle_diameter;
	// New re-windowed particle diameter (in Angstroms)
	DOUBLE new_particle_diameter;

	long int particle_radius2;

	// Low pass filetr cutoff for references (in Angstroms)
	DOUBLE lowpass;

	// Low pass filetr cutoff for particles (in Angstroms)
	DOUBLE lowpass_particle;

	// Original size of the reference images
	int particle_size;

	// Dimension of the filtered image
	int current_size;

	// Vector with all original reference images
	std::vector<MultidimArray<DOUBLE> > Mrefs;

	// FTs of the reference images for feature calculation
	std::vector<Projector > PPref;

	// Is density in micrograph inverted wrt templates?
	bool do_invert;

	// Modified by ZhouQ
	// Output ref-subtracted particles?
	bool do_outputrefsubptcls;
	// Root name of output ref-subtracted particles
	FileName fn_refsubptcls ;
	// Fit particles with references by least square method?
	// Does not work for small particles
	bool do_fitleastsquare ;
	// Fit particles scales?
	bool do_fitscale;

        // Perform normalization of the extract images
	bool do_normalise;

	// Subtract ramp instead of a level background in normalization
	bool do_ramp;

	// Standard deviations to remove black and white dust
	DOUBLE white_dust_stddev, black_dust_stddev;

	// Rewindow ref-subtracted particles?
	bool do_rewindowrefsub;

	// Rewindow raw particles?
	bool do_rewindowraw;

	// Root name of output rewindowed particles
	FileName fn_winptcls;
	
	// Only output star file containing updated location information of segmented particles?
	bool do_onlymakestar;
	
	// Remain coordinates of raw particles when only make star
	bool remain_rawcoord;
	
	// Update offset to segmented particles when only make star
	bool updateOffset2SegPart;
	// Root name of particles when only output star file
	FileName fn_onlystar;
	
	// Shift particles in Fourier transform?
	bool do_shiftinfouriertransform;

	// Tabulated sin and cosine functions for shifts in Fourier space
	TabSine tab_sin;
	TabCosine tab_cos;

	// The present symmetry from which expand the euler angles distribution to c1 symmetry
	FileName fn_symin;

	// Randomize present symmetry to C1 symmetry?
	bool do_randomsym;

	// The index of symmetry operation which is to be applied to segmented particles
	int index_sym;

	/** List of symmetry operators */
	std::vector <Matrix2D<DOUBLE> > R_repository, L_repository;

	// MultidimArray to store origin shift of a part of particle 
	MultidimArray<DOUBLE> rewinOriginOffsets;

	// Size of rewindowed particles
	int rewindow_size;
	// Coordinates of target center in pixel. The origin is the center of reference
	DOUBLE t_x, t_y, t_z;

	// Local Euler angles of target in original 3D map
	DOUBLE t_rot, t_tilt, t_psi;

	// Correct the references for CTF effects?
	bool do_ctf;

	// Keep the CTFs unchanged until the first peak?
	bool intact_ctf_first_peak;

public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// General function to decide what to do
	void run();

	// Initialise some general stuff after reading
	void initialise();

	void segmentOneParticle(long int ipart);
	// Modified by ZhouQ
	// fit Particle scale by least square method
	void fitParticle(Image<DOUBLE> &img, MultidimArray<DOUBLE> &Mref_rot);

protected:

	// Write out
	void write();

};


#endif /* PARTICLE_SEGMENTER_H_ */
