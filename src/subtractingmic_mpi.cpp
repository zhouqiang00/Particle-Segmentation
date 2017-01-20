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
#include "src/subtractingmic_mpi.h"

void SubtractingMicrographMpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    SubtractingMicrograph::read(argc, argv, node->rank);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? 1 : 0;

    // Possibly also read parallelisation-dependent variables here

    // Print out MPI info
	printMpiNodesMachineNames(*node);


}

void SubtractingMicrographMpi::initialise()
{
	// First read in non-parallelisation-dependent variables
	SubtractingMicrograph::initialise(node->rank,node->size);

}



void SubtractingMicrographMpi::run()
{

	// Each node does part of the work
	long int my_first_mic, my_last_mic, my_nr_mics;
	divide_equally(fn_allmics.size(), node->size, node->rank, my_first_mic, my_last_mic);
	my_nr_mics = my_last_mic - my_first_mic + 1;
	int barstep;
	if (do_renumberRawImages)
	{	
		if (verb > 0)
		{
			std::cout << "Renumbering all particles in the input star file..." << std::endl;
			init_progress_bar(my_nr_mics);
			barstep = XMIPP_MAX(1, my_nr_mics / 60);
		}
		for (int imic = 0; imic < fn_mics.size();imic++)
		{
			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
			renumberInputStarFileAndOutputRawImages(imic);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (node->isMaster())
		{
			joinStarFiles(fn_allmics, "_" + fn_renumber + ".star",fn_in.insertBeforeExtension("_" + fn_renumber));
		}
	}
	if (do_rewindowMicrographs)
	{
		if (verb > 0)
		{
			std::cout << "Processing all micrographs in the input star file..." << std::endl;
			init_progress_bar(my_nr_mics);
			barstep = XMIPP_MAX(1, my_nr_mics / 60);
		}
		for (int imic = 0; imic < fn_mics.size();imic++)
		{
			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
			subtractProjectionsFromMicrograph(imic);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if ((do_reextract || do_subtract_neighbor ) && node->isMaster())
		{
			joinStarFiles(fn_allmics, "_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_" + fn_rewindow));
		}
	}
	if (do_rewindowParticles)
	{
		if (verb > 0)
		{
			std::cout << "Processing all particles in the input star file..." << std::endl;
			init_progress_bar(my_nr_mics);
			barstep = XMIPP_MAX(1, my_nr_mics / 60);
		}
		for (int imic = 0; imic < fn_mics.size();imic++)
		{
			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
			subtractProjectionsFromParticles(imic);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (node->isMaster())
		{
			if (do_output_raw_particles)
			{
				joinStarFiles(fn_allmics, "_raw_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_raw_" + fn_rewindow) );
			}
			if (do_output_sub_particles)
			{
				joinStarFiles(fn_allmics, "_sub_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_sub_" + fn_rewindow) );
			}
			joinStarFiles(fn_allmics, "_coord_" + fn_rewindow + ".star",fn_in.insertBeforeExtension("_coord_" + fn_rewindow) );
		}
	}
	if (verb > 0)
		std::cout << " Done!" <<std::endl;

}


