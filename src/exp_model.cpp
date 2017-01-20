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
#include "src/exp_model.h"

void ExpOriginalParticle::addParticle(long int _particle_id, int _random_subset, int _order)
{
	// Keep random_subsets equal in each original particle
	if (random_subset != _random_subset)
		REPORT_ERROR("ExpOriginalParticle:addParticle: incompatible random subsets between particle and its original particle");
	particles_id.push_back(_particle_id);
	particles_order.push_back(_order);
}

long int Experiment::numberOfParticles(int random_subset)
{
	if (random_subset == 0)
		return particles.size();
	else
	{
		long int result = 0;
		for (long int i = 0; i < ori_particles.size(); i++)
		{
			if (ori_particles[i].random_subset == random_subset)
			{
				result += ori_particles[i].particles_id.size();
			}
		}
		return result;
	}
}

long int Experiment::numberOfOriginalParticles(int random_subset)
{
	if (random_subset == 0)
		return ori_particles.size();
	else if (random_subset == 1)
		return nr_ori_particles_subset1;
	else if (random_subset == 2)
		return nr_ori_particles_subset2;
	else
		REPORT_ERROR("ERROR: Experiment::numberOfOriginalParticles invalid random_subset: " + integerToString(random_subset));
}


long int Experiment::numberOfMicrographs()
{
	return micrographs.size();
}

long int Experiment::numberOfGroups()
{
	return groups.size();
}

long int Experiment::getMicrographId(long int part_id)
{
#ifdef DEBUG_CHECKSIZES
	if (part_id >= particles.size())
	{
		std::cerr<< "part_id= "<<part_id<<" particles.size()= "<< particles.size() <<std::endl;
		REPORT_ERROR("part_id >= particles.size()");
	}
#endif
	return (particles[part_id]).micrograph_id;
}

long int Experiment::getGroupId(long int part_id)
{
#ifdef DEBUG_CHECKSIZES
	if (part_id >= particles.size())
	{
		std::cerr<< "part_id= "<<part_id<<" particles.size()= "<< particles.size() <<std::endl;
		REPORT_ERROR("part_id >= particles.size()");
	}
#endif
	return (particles[part_id]).group_id;
}

int Experiment::getRandomSubset(long int part_id)
{
	return particles[part_id].random_subset;
}

MetaDataTable Experiment::getMetaDataImage(long int part_id)
{
	MetaDataTable result;
	result.addObject(MDimg.getObject(part_id));
	return result;
}

long int Experiment::addParticle(long int group_id,
		long int micrograph_id, int random_subset)
{

	if (group_id >= groups.size())
		REPORT_ERROR("Experiment::addImage: group_id out of range");

	if (micrograph_id >= micrographs.size())
		REPORT_ERROR("Experiment::addImage: micrograph_id out of range");

	ExpParticle particle;
	particle.id = particles.size();
	particle.group_id = group_id;
	particle.micrograph_id = micrograph_id;
	particle.random_subset = random_subset;
	// Push back this particle in the particles vector
	particles.push_back(particle);
	(micrographs[micrograph_id].particle_ids).push_back(particle.id);

	// Return the id in the particles vector
	return particle.id;

}

long int Experiment::addOriginalParticle(std::string part_name, int _random_subset)
{

	ExpOriginalParticle ori_particle;
	ori_particle.random_subset = _random_subset;
	ori_particle.name = part_name;
	long int id = ori_particles.size();
	ori_particles.push_back(ori_particle);

	// Return the id in the ori_particles vector
	return id;

}

long int Experiment::addGroup(std::string group_name)
{
	// Add new group to this Experiment
	ExpGroup group;
	group.id = groups.size(); // start counting groups at 0!
	group.name = group_name;

	// Push back this micrograph
	groups.push_back(group);

	// Return the id in the micrographs vector
	return group.id;

}

long int Experiment::addAverageMicrograph(std::string avg_mic_name)
{
	// Add new averagemicrograph to this Experiment
	AverageMicrograph micrograph;
	micrograph.id = average_micrographs.size();
	micrograph.name = avg_mic_name;

	// Push back this micrograph
	average_micrographs.push_back(micrograph);

	// Return the id in the micrographs vector
	return micrograph.id;

}

long int Experiment::addMicrograph(std::string mic_name)
{
	// Add new micrograph to this Experiment
	ExpMicrograph micrograph;
	micrograph.id = micrographs.size();
	micrograph.name = mic_name;

	// Push back this micrograph
	micrographs.push_back(micrograph);

	// Return the id in the micrographs vector
	return micrograph.id;

}

void Experiment::divideOriginalParticlesInRandomHalves(int seed)
{

	// Only do this if the random_subset of all original_particles is zero
	bool all_are_zero = true;
	bool some_are_zero = false;
	nr_ori_particles_subset1 = 0;
	nr_ori_particles_subset2 = 0;
	for (long int i = 0; i < ori_particles.size(); i++)
	{
		int random_subset = ori_particles[i].random_subset;
		if (random_subset != 0)
		{
			all_are_zero = false;
			// Keep track of how many particles there are in each subset
			if (random_subset == 1)
				nr_ori_particles_subset1++;
			else if (random_subset == 2)
				nr_ori_particles_subset2++;
			else
				REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
		}
		else
			some_are_zero = true;

		if (!all_are_zero && some_are_zero)
			REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: some random subset values are zero and others are not. They should all be zero, or all bigger than zero!");
	}

	//std::cerr << " all_are_zero= " << all_are_zero << " some_are_zero= " << some_are_zero << std::endl;

	if (all_are_zero)
	{
		// Only randomise them if the random_subset values were not read in from the STAR file
		srand(seed);
		for (long int i = 0; i < ori_particles.size(); i++)
		{
			int random_subset = rand() % 2 + 1;
			ori_particles[i].random_subset = random_subset; // randomly 1 or 2
			if (random_subset == 1)
				nr_ori_particles_subset1++;
			else if (random_subset == 2)
				nr_ori_particles_subset2++;
			else
				REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));

			// Loop over all particles in each ori_particle and set their random_subset
			for (long int j = 0; j < ori_particles[i].particles_id.size(); j++)
			{
				long int part_id = (ori_particles[i]).particles_id[j];
				{
					particles[part_id].random_subset = random_subset;
					MDimg.setValue(EMDL_PARTICLE_RANDOM_SUBSET, random_subset, part_id);
				}
			}
		}
	}
}

void Experiment::randomiseOriginalParticlesOrder(int seed, bool do_split_random_halves)
{
	//This static flag is for only randomize once
	static bool randomised = false;
	if (!randomised)
	{

		srand(seed);
		std::vector<ExpOriginalParticle> new_ori_particles;

		if (do_split_random_halves)
		{
			std::vector<long int> ori_particle_list1, ori_particle_list2;
			ori_particle_list1.clear();
			ori_particle_list2.clear();
			// Fill the two particle lists
			for (long int i = 0; i < ori_particles.size(); i++)
			{
				int random_subset = ori_particles[i].random_subset;
				if (random_subset == 1)
					ori_particle_list1.push_back(i);
				else if (random_subset == 2)
					ori_particle_list2.push_back(i);
				else
					REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
			}

			// Just a silly check for the sizes of the ori_particle_lists (to be sure)
			if (ori_particle_list1.size() != nr_ori_particles_subset1)
				REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid ori_particle_list1 size:" + integerToString(ori_particle_list1.size()) + " != " + integerToString(nr_ori_particles_subset1));
			if (ori_particle_list2.size() != nr_ori_particles_subset2)
				REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid ori_particle_list2 size:" + integerToString(ori_particle_list2.size()) + " != " + integerToString(nr_ori_particles_subset2));

			// Randomise the two particle lists
			std::random_shuffle(ori_particle_list1.begin(), ori_particle_list1.end());
			std::random_shuffle(ori_particle_list2.begin(), ori_particle_list2.end());

			// First fill new_ori_particles with the first subset, then with the second
			for (long int i = 0; i < ori_particle_list1.size(); i++)
				new_ori_particles.push_back(ori_particles[ori_particle_list1[i]]);
			for (long int i = 0; i < ori_particle_list2.size(); i++)
				new_ori_particles.push_back(ori_particles[ori_particle_list2[i]]);

		}
		else
		{

			// First fill in order
			std::vector<long int> ori_particle_list;
			ori_particle_list.resize(ori_particles.size());
			for (long int i = 0; i < ori_particle_list.size(); i++)
				ori_particle_list[i] = i;

			// Randomise
			std::random_shuffle(ori_particle_list.begin(), ori_particle_list.end());

			// Refill new_ori_particles
			for (long int i = 0; i < ori_particle_list.size(); i++)
				new_ori_particles.push_back(ori_particles[ori_particle_list[i]]);
		}

		ori_particles=new_ori_particles;
		randomised = true;

	}
}


void Experiment::expandToMovieFrames(FileName fn_data_movie, int verb)
{

//#define DEBUG_EXPAND
#ifdef DEBUG_EXPAND
	Timer timer;
	int tall = timer.setNew("ALL");
	int tread = timer.setNew("read");
	int tsort = timer.setNew("sort");
	int tmakevec = timer.setNew("makevec");
	int tselect = timer.setNew("select from MDimg");
	int tsearchmic = timer.setNew("seacrh mic");
	int tsearchgroup = timer.setNew("seacrh group");
	int tsearchori = timer.setNew("seacrh oriparticle");
	int taddpart = timer.setNew("add particle");
	int torderori = timer.setNew("order particles in oriparticle");
	timer.tic(tall);
#endif

	MetaDataTable MDmovie;
#ifdef DEBUG_EXPAND
	timer.tic(tread);
#endif
	MDmovie.read(fn_data_movie);

#ifdef DEBUG_EXPAND
	std::cerr << "now read in MDmovie" << std::endl;
	timer.toc(tread);
#endif
	if (!MDmovie.containsLabel(EMDL_MICROGRAPH_NAME) || !MDmovie.containsLabel(EMDL_PARTICLE_ORI_NAME))
		REPORT_ERROR("Experiment::expandToMovieFrames Error: movie metadata file does not contain rlnMicrographName as well as rlnOriginalParticleName");

	// Sort movie STAR file on EMDL_MICROGRAPH_NAME
#ifdef DEBUG_EXPAND
	timer.tic(tsort);
#endif
	MDmovie.newSort(EMDL_MICROGRAPH_NAME, false, true); // false=no reverse, true= do sort only on string after "@"
#ifdef DEBUG_EXPAND
	//MDmovie.write("sorted_movie.star");
	//MDimg.write("sorted_MDimg.star");
	//std::cerr << "Written sorted_movie.star and sorted_MDimg.star" << std::endl;
        std::cerr << "now sorted MDmovie" << std::endl;
	timer.toc(tsort);
#endif

	// Re-build new Experiment Exp_movie from scratch
	Experiment Exp_movie;

#ifdef DEBUG_EXPAND
	timer.tic(tmakevec);
#endif
	// Make a temporary vector of all image names in the current Experiment to gain speed
	std::vector<FileName> fn_curr_imgs, fn_curr_groups;
	std::vector<int> count_frames;
	std::vector<long int> pointer_current_idx;
	FileName fn_curr_img;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
	{
		MDimg.getValue(EMDL_IMAGE_NAME, fn_curr_img);
		long int group_id;
		MDimg.getValue(EMDL_MLMODEL_GROUP_NO, group_id);
		fn_curr_imgs.push_back(fn_curr_img);
		fn_curr_groups.push_back(groups[group_id-1].name);
		count_frames.push_back(0);
	}
#ifdef DEBUG_EXPAND
	timer.toc(tmakevec);
#endif

	if (verb > 0)
		init_progress_bar(MDmovie.numberOfObjects());
	int bar_step = MDmovie.numberOfObjects() / 60;

	FileName last_found_fn_mic="";
	long int last_found_mic_idx = -1, last2_found_mic_idx = -1;
	long i_object = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmovie)
	{
		long int group_id, mic_id, part_id;
		int my_random_subset, my_class;
		DOUBLE rot, tilt, psi, xoff, yoff;
		FileName fn_curr_img, group_name, fn_ori_part, fn_mic;
		MDmovie.getValue(EMDL_PARTICLE_ORI_NAME, fn_ori_part);

		bool have_found = false;
		// Find this particle in the current Experiment
		// Start searching from the beginning of the last found micrograph (because both MDimg and MDmovie are sorted on that)
		// Particles are not necessarily ordered in the micrographs, so therefore start searching at first idx of each new micrograph!

#ifdef DEBUG_EXPAND
		timer.tic(tselect);
#endif
		// Not starting at the beginning every time makes this scale ~linearly, instead of squared! :-)
		// rlnMicrographName needs to be sorted in both MDimg and MDmovie for this to work (which is now the case)
		//for (long int idx = 0; idx < fn_curr_imgs.size(); idx++)
		for (long int idx = last2_found_mic_idx + 1; idx < fn_curr_imgs.size(); idx++)
		{
			// Find a match
			if (fn_curr_imgs[idx] == fn_ori_part)
			{

				// Now get the angles from the current Experiment
				MDimg.getValue(EMDL_ORIENT_ROT, rot, idx);
				MDimg.getValue(EMDL_ORIENT_TILT, tilt, idx);
				MDimg.getValue(EMDL_ORIENT_PSI, psi, idx);
				MDimg.getValue(EMDL_ORIENT_ORIGIN_X, xoff, idx);
				MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, yoff, idx);
				MDimg.getValue(EMDL_PARTICLE_CLASS, my_class, idx);
				// Also get the random subset (if present)
				if (!MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, my_random_subset, idx))
					my_random_subset = 0;
				MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic, idx);
				if (last_found_fn_mic != fn_mic)
				{
					last_found_fn_mic = fn_mic;
					last2_found_mic_idx = last_found_mic_idx;
					last_found_mic_idx = idx;
				}

				// count how many frames are measured for each particle
				count_frames[idx]++;
				// Also keep track to which particle each image in MDmovie belongs
				pointer_current_idx.push_back(idx);
				group_name = fn_curr_groups[idx];
				have_found = true;
				break;
			}

		}
#ifdef DEBUG_EXPAND
		timer.toc(tselect);
#endif

		// Only include particles that were already in the current Experiment
		if (have_found)
		{
			// Add new micrographs or get mic_id for existing micrograph
			FileName mic_name;
			MDmovie.getValue(EMDL_MICROGRAPH_NAME, mic_name);

			// If this micrograph did not exist in the Exp_movie yet, add it to the Exp_movie experiment
#ifdef DEBUG_EXPAND
		timer.tic(tsearchmic);
#endif
			mic_id = -1;
			for (long int i = Exp_movie.micrographs.size() - 1; i >= 0; i--)
			{
				if (Exp_movie.micrographs[i].name == mic_name)
				{
					mic_id = Exp_movie.micrographs[i].id;
					break;
				}
			}
			if (mic_id < 0)
				mic_id = Exp_movie.addMicrograph(mic_name);
#ifdef DEBUG_EXPAND
		timer.toc(tsearchmic);
#endif

			// Add frameno@ to existing group names, so that separate weighting may be applied to different dose images
			// NO THIS HAS NO SENSE IF WE'RE ONLY DOING ONE ITERATION ANYWAY!!! THEN IT'S JUST A WASTE OF MEMORY....

#ifdef DEBUG_EXPAND
		timer.tic(tsearchgroup);
#endif
			// If this group did not exist yet, add it to the experiment
			group_id = -1;
			for (long int i = Exp_movie.groups.size() - 1; i >= 0; i--)
			{
				if (Exp_movie.groups[i].name == group_name)
				{
					group_id = Exp_movie.groups[i].id;
					break;
				}
			}
			if (group_id < 0)
				group_id = Exp_movie.addGroup(group_name);
#ifdef DEBUG_EXPAND
		timer.toc(tsearchgroup);
#endif

#ifdef DEBUG_EXPAND
		timer.tic(taddpart);
#endif
			// Create a new particle
			part_id = Exp_movie.addParticle(group_id, mic_id, my_random_subset);
                        Exp_movie.MDimg.addObject();

			// Copy the current row of MDimgin into the current row of MDimg
			Exp_movie.MDimg.setObject(MDmovie.getObject(), part_id);

			// Set the orientations
			Exp_movie.MDimg.setValue(EMDL_ORIENT_ROT, rot, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_TILT, tilt, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_PSI, psi, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_ORIGIN_X, xoff, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y, yoff, part_id);
			// Now also set the priors on the orientations equal to the orientations from the averages
			Exp_movie.MDimg.setValue(EMDL_ORIENT_ROT_PRIOR, rot, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_TILT_PRIOR, tilt, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_PSI_PRIOR, psi, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_ORIGIN_X_PRIOR, xoff, part_id);
			Exp_movie.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, yoff, part_id);
			Exp_movie.MDimg.setValue(EMDL_PARTICLE_CLASS, my_class, part_id);
			Exp_movie.MDimg.setValue(EMDL_PARTICLE_RANDOM_SUBSET, my_random_subset, part_id);
			// Set normcorrection to 1
			DOUBLE norm = 1.;
			Exp_movie.MDimg.setValue(EMDL_IMAGE_NORM_CORRECTION, norm);
			// Set the new group name
			Exp_movie.MDimg.setValue(EMDL_MLMODEL_GROUP_NAME, group_name);
#ifdef DEBUG_EXPAND
		timer.toc(taddpart);
#endif

#ifdef DEBUG_EXPAND
		timer.tic(tsearchori);
#endif
			// Add ExpOriParticles
			// If this ori_particle did not exist in the Exp_movie yet, add it to the Exp_movie experiment
			long int ori_part_id = -1;
			// Loop backward to find the previous ori_particle first
			for (long int i = Exp_movie.ori_particles.size() - 1; i >= 0; i--)
			{
				if (Exp_movie.ori_particles[i].name == fn_ori_part)
				{
					ori_part_id = i;
					break;
				}
			}
			// If no ExpOriParticles with this name was found, then add new one
			if (ori_part_id < 0)
				ori_part_id = Exp_movie.addOriginalParticle(fn_ori_part, my_random_subset);
#ifdef DEBUG_EXPAND
		timer.toc(tsearchori);
#endif


			// Add this particle to the OriginalParticle
			// get Number from mic_name (-1 if empty mic_name, or no @ in mic_name)
			std::string fnt;
			long int my_order;
			mic_name.decompose(my_order, fnt);
			(Exp_movie.ori_particles[ori_part_id]).addParticle(part_id, my_random_subset, my_order);

		}

		if (verb > 0 && i_object % bar_step == 0)
			progress_bar(i_object);
		i_object++;
	}

	if (verb > 0)
		progress_bar(MDmovie.numberOfObjects());

	if (Exp_movie.MDimg.numberOfObjects() == 0)
		REPORT_ERROR("Experiment::expandToMovieFrames: ERROR: no movie frames selected. Check filenames of micrographs, movies and particle stacks!");

	// Now that all particles from MDmovie have been parsed, set nr_frames per particle in the metadatatable
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(Exp_movie.MDimg)
	{
		Exp_movie.MDimg.setValue(EMDL_PARTICLE_NR_FRAMES, count_frames[(pointer_current_idx[current_object])]);
	}

	// Now replace the current Experiment with Exp_movie
	(*this) = Exp_movie;

#ifdef DEBUG_EXPAND
		timer.tic(torderori);
#endif
	// Order the particles in each ori_particle
	orderParticlesInOriginalParticles();
#ifdef DEBUG_EXPAND
		timer.toc(torderori);
#endif

#ifdef DEBUG_EXPAND
	timer.toc(tall);
	timer.printTimes(true);
#endif


}
void Experiment::orderParticlesInOriginalParticles()
{
	// If the orders are negative (-1) then dont sort anything
	if (ori_particles[0].particles_order[0] < 0)
		return;

	for (long int i = 0; i < ori_particles.size(); i++)
	{
		int nframe = ori_particles[i].particles_order.size();

		std::vector<std::pair<long int, long int> > vp;
        vp.reserve(nframe);
        for (long int j = 0; j < nframe; j++)
        	vp.push_back(std::make_pair(ori_particles[i].particles_order[j], j));
        // Sort on the first elements of the pairs
        std::sort(vp.begin(), vp.end());

        // tmp copy of particles_id
        std::vector<long int> _particles_id = ori_particles[i].particles_id;
        for (int j = 0; j < nframe; j++)
			ori_particles[i].particles_id[j] = _particles_id[vp[j].second];

		// We now no longer need the particles_order vector, clear it to save memory
		ori_particles[i].particles_order.clear();
	}

}

void Experiment::usage()
{
	std::cout
	<< "  -i                     : Starfile with input images\n"
	;
}

// Read from file
void Experiment::read(FileName fn_exp, bool do_ignore_original_particle_name, bool do_ignore_group_name, bool do_preread_images)
{

//#define DEBUG_READ
#ifdef DEBUG_READ
	std::cerr << "Entering Experiment::read" << std::endl;
	Timer timer;
	int tall = timer.setNew("ALL");
	int tread = timer.setNew("read");
	int tsort = timer.setNew("sort");
	int tfill = timer.setNew("fill");
	int tgroup = timer.setNew("find group");
	int tori = timer.setNew("find ori_particle");
	int tdef = timer.setNew("set defaults");
	int tend = timer.setNew("ending");
	char c;
	timer.tic(tall);
	timer.tic(tread);
#endif

	// Initialize by emptying everything
	clear();
	long int group_id, mic_id, part_id;

	if (!fn_exp.isStarFile())
	{
		// Read images from stack. Ignore all metadata, just use filenames
		// Add a single Micrograph
		group_id = addGroup("group");
		mic_id = addMicrograph("micrograph");

		// Check that a MRC stack ends in .mrcs, not .mrc (which will be read as a MRC 3D map!)
		if (fn_exp.contains(".mrc") && !fn_exp.contains(".mrcs"))
			REPORT_ERROR("Experiment::read: ERROR: MRC stacks of 2D images should be have extension .mrcs, not .mrc!");

		// Read in header-only information to get the NSIZE of the stack
		Image<DOUBLE> img;
		img.read(fn_exp, false); // false means skip data, only read header
		
		// allocate 1 block of memory
		particles.reserve(NSIZE(img()));
		ori_particles.reserve(NSIZE(img()));
		for (long int n = 0; n <  NSIZE(img()); n++)
		{
			FileName fn_img;
			fn_img.compose(n+1, fn_exp); // fn_img = integerToString(n) + "@" + fn_exp;
			// Add the particle to my_area = 0
			part_id = addParticle(group_id, mic_id);
                        MDimg.addObject();
			if (do_preread_images)
			{
				Image<DOUBLE> img;
				img.read(fn_img);
				img().setXmippOrigin();
				particles[part_id].img = img();
			}
			// Also add OriginalParticle
			(ori_particles[addOriginalParticle("particle")]).addParticle(part_id, 0, -1);
			// Set the filename and other metadata parameters
			MDimg.setValue(EMDL_IMAGE_NAME, fn_img, part_id);
		}

	}
	else
	{
		// Just read first data block
		MDimg.read(fn_exp);

#ifdef DEBUG_READ
	std::cerr << "Done reading MDimg" << std::endl;
	timer.toc(tread);
	timer.tic(tsort);
	//std::cerr << "Press any key to continue..." << std::endl;
	//std::cin >> c;
#endif

		// Sort input particles on micrographname
		bool is_mic_a_movie = false, star_contains_micname;
		star_contains_micname = MDimg.containsLabel(EMDL_MICROGRAPH_NAME);
		if (star_contains_micname)
		{
			// See if the micrograph names contain an "@", i.e. whether they are movies and we are inside polishing or so.
			FileName fn_mic;
			MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			if (fn_mic.contains("@"))
			{
				is_mic_a_movie = true;
				MDimg.newSort(EMDL_MICROGRAPH_NAME, false, true); // sort on part AFTER "@"
			}
			else
			{
				is_mic_a_movie = false;
				MDimg.newSort(EMDL_MICROGRAPH_NAME); // just sort on fn_mic
			}

			if (do_ignore_group_name)
				group_id = addGroup("group");
		}
		else
		{
			// If there is no EMDL_MICROGRAPH_NAME, then just use a single group and micrograph
			group_id = addGroup("group");
			mic_id = addMicrograph("micrograph");
		}
#ifdef DEBUG_READ
	std::cerr << "Done sorting MDimg" << std::endl;
	timer.toc(tsort);
	timer.tic(tfill);
	long nr_read = 0;
#endif
                // allocate 1 block of memory
                particles.reserve(MDimg.numberOfObjects());
               
		// Now Loop over all objects in the metadata file and fill the logical tree of the experiment
		long int last_oripart_idx = -1;
		int nr_frames = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
		{
			// Add new micrographs or get mic_id for existing micrograph
  	                mic_id = -1;
                        long int idx = micrographs.size();
			long int avg_mic_idx = average_micrographs.size() - 1;
                        std::string last_mic_name = (idx > 0) ? micrographs[idx-1].name : "";
			FileName mic_name=""; // Filename instead of string because will decompose below
			std::string mic_name_after_at="", group_name="", last_mic_name_after_at="", last_ori_mic_name="";
			if (is_mic_a_movie)
				last_mic_name_after_at = (idx > 0) ? last_mic_name.substr(last_mic_name.find("@")+1) : "";			

			if (star_contains_micname)
			{
				MDimg.getValue(EMDL_MICROGRAPH_NAME, mic_name);

				// Now the micrograph names are sorted, but for movies only the part after "@" is sorted....
				long idx = micrographs.size();
				if (is_mic_a_movie)
				{
					// Go backward in old micrographs and see if this frame from the movie has been seen already
					// But only for the same mic_name_after_at!
					mic_name_after_at = mic_name.substr(mic_name.find("@")+1);
					if (last_mic_name_after_at == mic_name_after_at)
					{
						// TODO: only go back maximum of nr_frames steps!
						for (long i = micrographs.size() - 1; i >= 0; i--)
						{
							// if we find this frame from this mic_name, then set mic_id and break
							if (micrographs[i].name == mic_name)
							{
								mic_id = micrographs[i].id;
								break;
							}
						}
					}
					else
					{
						// A new (original) micrograph
						last_oripart_idx = ori_particles.size();
						avg_mic_idx = -1;
						// How many movie frames per micrograph are there? (Only count for first average_micrograph, all the rest is the same)
						if (average_micrographs.size() == 1)
						{
							nr_frames = micrographs.size();
							ori_particles.reserve(MDimg.numberOfObjects()/nr_frames);
#ifdef DEBUG_READ							
							std::cerr << "nr_frames = " << nr_frames << std::endl;
#endif
						}
					}
				}
				else // not a movie
				{
					if (last_mic_name == mic_name)
					{
						// This particle belongs to the previous micrograph
						mic_id = micrographs[idx - 1].id;
					}
					else
					{
						// A new micrograph
						last_oripart_idx = ori_particles.size();
					}
				}

				// Make a new micrograph
				if (mic_id < 0)
					mic_id = addMicrograph(mic_name);

				// Make a new average_micrograph
				if (is_mic_a_movie && avg_mic_idx < 0)
				{
					avg_mic_idx = addAverageMicrograph(mic_name_after_at);
				}

#ifdef DEBUG_READ
				timer.tic(tgroup);
#endif

				// For example in particle_polishing the groups are not needed...
				if (!do_ignore_group_name)
				{
					// Check whether there is a group label, if not use a group for each micrograph
					if (MDimg.containsLabel(EMDL_MLMODEL_GROUP_NAME))
						MDimg.getValue(EMDL_MLMODEL_GROUP_NAME, group_name);
					else
						group_name = mic_name;

					// If this group did not exist yet, add it to the experiment
					group_id = -1;
					for (long int i = groups.size() - 1; i >= 0; i--) // search backwards to find match faster
					{
						if (groups[i].name == group_name)
						{
							group_id = groups[i].id;
							break;
						}
					}
					if (group_id < 0)
						group_id = addGroup(group_name);
				}

#ifdef DEBUG_READ
				timer.toc(tgroup);
#endif

			}
			else
			{
				// All images belong to the same micrograph
				mic_id = 0;
				group_id = 0;
			}

			// If there is an EMDL_PARTICLE_RANDOM_SUBSET entry in the input STAR-file, then set the random_subset, otherwise use default (0)
			int my_random_subset;
			if (!MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, my_random_subset))
				my_random_subset = 0;

			// Create a new particle
			part_id = addParticle(group_id, mic_id, my_random_subset);

#ifdef DEBUG_READ
			timer.tic(tori);
#endif

                        if (do_preread_images)
                        {
                            FileName fn_img;
                            MDimg.getValue(EMDL_IMAGE_NAME, fn_img);
                            Image<DOUBLE> img;
                            img.read(fn_img);
                            img().setXmippOrigin();
                            particles[part_id].img = img();
                        }

			// Add this particle to an existing OriginalParticle, or create a new OriginalParticle
			std::string ori_part_name;
			long int ori_part_id = -1;

			if (MDimg.containsLabel(EMDL_PARTICLE_ORI_NAME))
				MDimg.getValue(EMDL_PARTICLE_ORI_NAME, ori_part_name);
			else
				MDimg.getValue(EMDL_IMAGE_NAME, ori_part_name);

			if (MDimg.containsLabel(EMDL_PARTICLE_ORI_NAME) && !do_ignore_original_particle_name)
			{
				// Only search ori_particles for the last (original) micrograph
				for (long int i = last_oripart_idx; i < ori_particles.size(); i++)
				{
					if (ori_particles[i].name == ori_part_name)
					{
						ori_part_id = i;
						break;
					}
				}
			}

			// If no OriginalParticles with this name was found,
			// or if no EMDL_PARTICLE_ORI_NAME in the input file, or if do_ignore_original_particle_name
			// then add a new ori_particle
			if (ori_part_id < 0)
			{
				ori_part_id = addOriginalParticle(ori_part_name, my_random_subset);
				// Also add this original_particle to an original_micrograph (only for movies)
				if (is_mic_a_movie)
				{	
					average_micrographs[avg_mic_idx].ori_particles_id.push_back(ori_part_id);
				}
			}
#ifdef DEBUG_READ
			timer.toc(tori);
#endif


			// Add this particle to the OriginalParticle
			std::string fnt;
			long int my_order;
			mic_name.decompose(my_order, fnt);
			(ori_particles[ori_part_id]).addParticle(part_id, my_random_subset, my_order);

			// The group number is only set upon reading: it is not read from the STAR file itself,
			// there the only thing that matters is the order of the micrograph_names
			// Write igroup+1, to start numbering at one instead of at zero
			MDimg.setValue(EMDL_MLMODEL_GROUP_NO, group_id + 1, part_id);

#ifdef DEBUG_READ
			nr_read++;
#endif
		} // end loop over all objects in MDimg

#ifdef DEBUG_READ
		timer.toc(tfill);
		timer.tic(tdef);
		std::cerr << " MDimg.lastObject()= " << MDimg.lastObject() << std::endl;
		std::cerr << " nr_read= " << nr_read << " particles.size()= " << particles.size() << " ori_particles.size()= " << ori_particles.size()  << " micrographs.size()= " << micrographs.size() << " average_micrographs.size()= " << average_micrographs.size() << " groups.size()= " << groups.size() << std::endl;
#endif

	}

#ifdef DEBUG_READ
	std::cerr << "Done filling MDimg" << std::endl;
	//std::cerr << "Press any key to continue..." << std::endl;
	//std::cin >> c;
#endif

	// Make sure some things are always set in the MDimg
	bool have_rot  = MDimg.containsLabel(EMDL_ORIENT_ROT);
	bool have_tilt = MDimg.containsLabel(EMDL_ORIENT_TILT);
	bool have_psi  = MDimg.containsLabel(EMDL_ORIENT_PSI);
	bool have_xoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_X);
	bool have_yoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Y);
	bool have_clas = MDimg.containsLabel(EMDL_PARTICLE_CLASS);
	bool have_norm = MDimg.containsLabel(EMDL_IMAGE_NORM_CORRECTION);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
	{
		DOUBLE dzero=0., done=1.;
		int izero = 0;
		if (!have_rot)
			MDimg.setValue(EMDL_ORIENT_ROT, dzero);
		if (!have_tilt)
			MDimg.setValue(EMDL_ORIENT_TILT, dzero);
		if (!have_psi)
			MDimg.setValue(EMDL_ORIENT_PSI, dzero);
		if (!have_xoff)
			MDimg.setValue(EMDL_ORIENT_ORIGIN_X, dzero);
		if (!have_yoff)
			MDimg.setValue(EMDL_ORIENT_ORIGIN_Y, dzero);
		if (!have_clas)
			MDimg.setValue(EMDL_PARTICLE_CLASS, izero);
		if (!have_norm)
			MDimg.setValue(EMDL_IMAGE_NORM_CORRECTION, done);
	}

#ifdef DEBUG_READ
	timer.toc(tdef);
	std::cerr << "Done setting defaults MDimg" << std::endl;
	timer.tic(tend);
	//std::cerr << "Press any key to continue..." << std::endl;
	//std::cin >> c;
#endif
	
	// Also set the image_size (use the last image for that, still in fn_img)
	FileName fn_img;
	Image<DOUBLE> img;
	MDimg.getValue(EMDL_IMAGE_NAME, fn_img, MDimg.firstObject());
	if (fn_img != "")
	{
		img.read(fn_img, false); //false means read only header, skip real data
		int image_size = XSIZE(img());
		if (image_size != YSIZE(img()))
			REPORT_ERROR("Experiment::read: xsize != ysize: only squared images allowed");
			// Add a single object to MDexp
			MDexp.addObject();
		MDexp.setValue(EMDL_IMAGE_SIZE, image_size);
		if (ZSIZE(img()) > 1)
		{
			if (image_size != ZSIZE(img()))
				REPORT_ERROR("Experiment::read: xsize != zsize: only cubed images allowed");
			MDexp.setValue(EMDL_IMAGE_DIMENSIONALITY, 3);
		}
		else
		{
			MDexp.setValue(EMDL_IMAGE_DIMENSIONALITY, 2);
		}
	}
	else
	{
		REPORT_ERROR("There are no images read in: please check your input file...");
	}

	// Order the particles in each ori_particle (only useful for realignment of movie frames)
	orderParticlesInOriginalParticles();

#ifdef DEBUG_READ
	timer.toc(tend);
	timer.toc(tall);
	timer.printTimes(false);
	//std::cerr << "Writing out debug_data.star" << std::endl;
	//write("debug");
	//exit(0);
#endif
}

// Write to file
void Experiment::write(FileName fn_root)
{

	std::ofstream  fh;
	FileName fn_tmp = fn_root+"_data.star";
    fh.open((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"Experiment::write: Cannot write file: " + fn_tmp);

    // Always write MDimg
    MDimg.write(fh);

	fh.close();

}
