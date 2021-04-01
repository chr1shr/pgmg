/** \file multigrid.cc
 * \brief Function implementations for the multigrid class. */

#include "multigrid.hh"

/** The multigrid destructor frees the dynamically allocated regions. */
multigrid::~multigrid() {
	for(int i=nlevels-1;i>=0;i--) {
		if(geo[i]!=NULL) delete geo[i];
		delete reg[i];
	}
	delete [] reg;
	delete [] geo;
}

/** Allocates the lower levels of a multigrid hierarchy. */
void multigrid::setup_hierarchy(geometry &gm) {
	int cproc=gm.procs;
	bool regular=false;

	// Set up the array for additional geometries. Since the top level will
	// never have a new geometry, assign it a null pointer.
	*(geo=new geometry*[mg_max_levels])=NULL;

	// Allocate the grid hierarchy, until the number of grid points
	while(reg[nlevels]->total()>1024) {

		if(nlevels==mg_max_levels)
		       p_fatal_error("Maximum number of multigrid levels reached",1);

		// If the minimum size is less than a tolerance, then reallocate the grid
		if(regular&&cproc>1&&reg[nlevels]->min_size()<20) {
			regular=false;

			// Compute a new child processor geometry, and tell the
			// parent regions to set up a table to communicate to
			// the new geometry
			geo[nlevels+1]=new geometry(reg[nlevels],cproc,8);
			cproc=geo[nlevels+1]->procs;
			reg[nlevels]->setup_outgoing_table(*(geo[nlevels+1]));

			// Check to see if this processor is not involved
			// in the child geometry
			if(geo[nlevels+1]->rank==-1) {

				// Delete the child geometry class, since it is
				// not used on this processor
				delete geo[nlevels+1];
				geo[nlevels+1]=NULL;

				// Check that there is enough memory for
				// sending the transfer. This processor only
				// needs to send data and doesn't need to
				// allocate space for receiving messages.
				reg[nlevels]->check_transfer_buf();

				// Send the RAT bounds
				reg[nlevels]->send_S();
				reg[nlevels]->send_wait();
				break;
			}
			reg[nlevels+1]=new region(reg[nlevels],*geo[nlevels+1],com);
		} else {

			// Allocate a new lower-level grid
			regular=true;
			reg[nlevels+1]=new region(reg[nlevels],com);
			geo[nlevels+1]=NULL;
		}
		nlevels++;
	}
	nlevels++;

	// Find the maximum number of levels
	MPI_Allreduce(&nlevels,&mlevels,1,MPI_INT,MPI_MAX,gm.cart);
#if PGMG_VERBOSE == 3
	printf("MG hierarchy: rank=%d, nlevels=%d, max levels=%d\n",rank,nlevels,mlevels);
#endif

	// Try and enable exact solution on the bottom level
	if(nlevels==mlevels) {
            reg[nlevels-1]->enable_exact();
            if(reg[nlevels-1]->Aexact!=NULL) reg[nlevels-1]->lu_decomposition();
        }
}

/** Performs a multigrid V-cycle. */
void multigrid::v_cycle(int mode) {
    int i,l;
    // Do downwards Gauss-Seidel operations and restrictions
    for(i=0;i<mlevels-1;i++) {
        for(l=0;l<v_down+i;l++) gauss_seidel(mode, i);
        restriction(i+1);
    }

    // Carry out Gauss-Seidel operations on the bottom level
    solve_exact(mlevels-1);

    // Do upwards Gauss-Seidel operations and interpolations
    for(i=mlevels-1;i>0;i--) {
        interpolation(i);
        for(l=0;l<v_up+i-1;l++) gauss_seidel(mode, i-1);
    }

    communicate();
}

/** Performs a full multigrid cycle. */
void multigrid::fmg_cycle() {
	int i,j,l;
	int mode=0;
	// Do downwards Gauss-Seidel operations and restrictions
	for(i=0;i<mlevels-1;i++) {
		restriction(i+1);
		for(l=0;l<v_down;l++) gauss_seidel(mode, i+1);
	}

	for(j=mlevels-1;j>0;j--) {
		for(i=j;i<mlevels-1;i++) {
			restriction(i+1);
			for(l=0;l<v_down;l++) gauss_seidel(mode, i+1);
		}

		// Carry out Gauss-Seidel operations on the bottom level
		for(l=0;l<v_bottom;l++) gauss_seidel(mode, mlevels-1);

		// Do upwards Gauss-Seidel operations and interpolations
		for(i=mlevels-1;i>=j;i--) {
			interpolation(i);
			for(l=0;l<v_up;l++) gauss_seidel(mode, i-1);
		}
	}
	communicate();
}

/** Computes the RAT matrices for the lower levels, once the top level matrices
 * have been initialized. */
void multigrid::compute_rats() {
	int i;
	for(i=1;i<nlevels;i++) {
		if(geo[i]==NULL) reg[i]->compute_rats();
		else {
			reg[i-1]->send_A();
			reg[i]->receive_A();
			reg[i-1]->send_wait();
		}
	}
	if(i<mlevels) {
		reg[i-1]->send_A();
		reg[i-1]->send_wait();
	}
}

/** Outputs a z-slice of the x field at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
void multigrid::output_x(const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_x(filename,k);
}

/** Outputs a z-slice of the r field at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
void multigrid::output_r(const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_r(filename,k);
}

/** Outputs a z-slice of the residual at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
void multigrid::output_residual(const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_residual(filename,k);
}

/** Carries out a Gauss-Seidel sweep at a certain level.
 * \param[in] level the region level to consider. */
void multigrid::gauss_seidel(int mode, int level) {
	if(level<nlevels&&reg[level]->gs_enable) reg[level]->gauss_seidel(mode);
}

/** Communicate the solution after cycles are done.
 * \param[in] level the region level to consider.
 * 			  Default is top level, 0.
 */
void multigrid::communicate(int level){
	if(level<nlevels) reg[level]->communicate();
}


/** Carries out an interpolation, interpolating the x field at the ith level
 * and adding it to the x field at the (i-1)th level.
 * \param[in] level the region level to consider. */
void multigrid::interpolation(int level) {
	if(level<nlevels) {
		if(geo[level]==NULL) reg[level]->interpolation();
		else {
			reg[level]->send_x();
			reg[level-1]->receive_x();
			reg[level]->send_wait();
		}
	} else if(level==nlevels)
		reg[level-1]->receive_x();
}

/** Carries out a restriction, computing the residual at the ith level and
 * passing it into the r field at the (i+1)th level.
 * \param[in] level the region level to consider. */
void multigrid::restriction(int level) {
	if(level<nlevels) {
		if(geo[level]==NULL) reg[level]->restriction();
		else {
			reg[level-1]->send_r();
			reg[level]->receive_r();
			reg[level-1]->send_wait();
		}
		reg[level]->clear_x_field();
	} else if(level==nlevels) {
		reg[level-1]->send_r();
		reg[level-1]->send_wait();
	}
}

/** Computes the max norm error at a given level. Allreduce, not just reduce.
 * \param[in] level the region level to consider.
 * \return The global error if this is the zeroth processor, the local error if
 * this processor is involved at the given level, and zero otherwise. */
double multigrid::l0_error_all(int level) {
	return nlevels>level?reg[level]->l0_error_all():0;
}

/** Computes the max norm error at a given level.
 * \param[in] level the region level to consider.
 * \return The global error if this is the zeroth processor, the local error if
 * this processor is involved at the given level, and zero otherwise. */
double multigrid::l0_error(int level) {
	return nlevels>level?reg[level]->l0_error():0;
}

/** Computes the L2 error at a given level. Allreduce, not just reduce.
 * \param[in] level the region level to consider.
 * \return The global error if this is the zeroth processor, the local error if
 * this processor is involved at the given level, and zero otherwise. */
double multigrid::l2_error_all(int level) {
	return nlevels>level?reg[level]->l2_error_all():0;
}

/** Computes the L2 error at a given level.
 * \param[in] level the region level to consider.
 * \return The global error if this is the zeroth processor, the local error if
 * this processor is involved at the given level, and zero otherwise. */
double multigrid::l2_error(int level) {
	return nlevels>level?reg[level]->l2_error():0;
}

/** Compute the L2 norm (not sqrt) of the right hand side */
double multigrid::rhs_l2norm(int level) {
	return nlevels>level?reg[level]->rhs_l2norm():0;
}

/**
 * update the predictor struct index, taking into
 * account the number of extra iterations needed to
 * reach tolerance above the previous timestep's prediction
 */
void multigrid::mg_guide::update(int extras) {

	// reduce iter value by 1 and increase to reflect any
	// extra iterations on previous step
	int diter = -1 + period*extras;

	// if using a moving average...
	if (use_avg) {

		// calculate any deviation from the moving average
		double di=iter-iter_avg;

		// reduce the index for each standard deviation by
		// which the average is exceeded
		int zscore = static_cast<int>(di/sqrt(iter_var));
		if (zscore > 0) diter -= zscore;

		// update averages
		iter_avg = (1-alpha)*iter_avg + alpha*iter;
		iter_var = (1-alpha)*(iter_var + alpha*di*di);
	}

	// update iteration index (reset if below period)
	iter += diter;
	if (iter < period) iter=2*period-1;
}
