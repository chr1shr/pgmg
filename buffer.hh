/** \file buffer.hh
 * \brief Header and function implementations for the communication buffer
 * class. */

#ifndef PGMG_BUFFER_HH
#define PGMG_BUFFER_HH

#include <cstdlib>
#include <cstdio>

/** \brief A class that manages common buffer space.
 *
 * Each processor has one of these classes that is commonly used by all regions
 * on that processor. The class has routines to scale up the size of the buffer
 * as necessary. */
class comm_buffer {
	public:
		/** A pointer the common buffer space. Since most of the buffer
		 * operations make use of doubles, keep it as a double pointer
		 * that can be cast to other types. */
		double *buf;
		/** The current size (in bytes) of the common buffer. */
		size_t mem;
		comm_buffer() : mem(0) {};
		~comm_buffer() {
			if(mem>0) free(buf);
		}
		inline void check_buf(int q) {reallocate(q*sizeof(double));}
		inline void check_buf_int(int q) {reallocate(q*sizeof(int));}
		inline void check_buf_uint(int q) {reallocate(q*sizeof(unsigned int));}
		inline void check_buf_float(int q) {reallocate(q*sizeof(float));}
		inline void check_buf_mesh(int q) {reallocate_mesh(q*sizeof(double));}
		inline void check_buf_int_mesh(int q) {reallocate_mesh(q*sizeof(int));}
		inline void check_buf_intstar_mesh(int q) {reallocate_mesh(q*sizeof(int*));}
		inline void check_buf_float_mesh(int q) {reallocate_mesh(q*sizeof(float));}
		inline void check_buf_mesh2(int q) {reallocate_mesh2(q*sizeof(double));}
		inline void check_buf_int_mesh2(int q) {reallocate_mesh2(q*sizeof(int));}
		inline void check_buf_intstar_mesh2(int q) {reallocate_mesh2(q*sizeof(int*));}
		inline void check_buf_float_mesh2(int q) {reallocate_mesh2(q*sizeof(float));}
	private:
		void reallocate(size_t nmem);
		void reallocate_mesh(size_t nmem);
		void reallocate_mesh2(size_t nmem);
};

#endif
