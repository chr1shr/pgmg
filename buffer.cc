#include "buffer.hh"
#include "common.hh"

void comm_buffer::reallocate(size_t nmem) {
	if(nmem>mem) {
		if(mem>0) free(buf);
		mem=nmem;
		buf=(double*) malloc(mem);
		if(buf==NULL) p_fatal_error("Buffer allocation failed",PGMG_MEMORY_ERROR);
	}
}
void comm_buffer::reallocate_mesh(size_t nmem) {
	//printf("Combuf %d\n",int(nmem));fflush(stdout);
	if(nmem>mem) {
		if(mem>0) free(buf);
                  while (mem < nmem) {
                    // Double memory
                    if (mem>0){
                      mem = 2*mem;
                    }
                    else mem=nmem;
                  }
                  if (mem<20480) mem=20480;
		buf=(double*) malloc(mem);
		if(buf==NULL) p_fatal_error("Buffer allocation failed",PGMG_MEMORY_ERROR);
	}
	else if(nmem<mem/2) {
		if(mem>0) free(buf);
                  if (nmem==0) mem=0;
                  else {
                    while (nmem < mem/2) {
                      // Halve memory
                      mem = mem/2;
                    }
                  }
                  if (mem<20480) mem=20480;
		buf=(double*) malloc(mem);
		if(buf==NULL) p_fatal_error("Buffer allocation failed",PGMG_MEMORY_ERROR);
	}
}
void comm_buffer::reallocate_mesh2(size_t nmem) {
	//printf("Combuf %d\n",int(nmem));fflush(stdout);
	if(nmem>mem) {
		if(mem>0) free(buf);
		  while (mem < nmem) {
		    // Double memory
		    if (mem>0){
		      mem = 2*mem;
		    }
		    else mem=nmem;
		  }
		  if (mem<20480) mem=20480;
		buf=(double*) malloc(mem);
		if(buf==NULL) p_fatal_error("Buffer allocation failed",PGMG_MEMORY_ERROR);
	}
}
