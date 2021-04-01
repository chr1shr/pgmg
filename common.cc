#include <cstdlib>

#include "common.hh"

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* safe_fopen(const char* filename,const char* mode) {
	FILE *temp=fopen(filename,mode);
	if(temp==NULL) {
		fprintf(stderr,"Error opening file \"%s\"",filename);
		exit(1);
	}
	return temp;
}

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* p_safe_fopen(const char* filename,const char* mode) {
	FILE *fp=fopen(filename,mode);
	if(fp==NULL) {
		fprintf(stderr,"Unable to open file '%s'\n",filename);
		MPI_Abort(world,PGMG_FILE_ERROR);
	}
	return fp;
}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void fatal_error(const char *p,int code) {
	fprintf(stderr,"Error: %s\n",p);
	exit(code);
}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void p_fatal_error(const char *p,int status) {
	fprintf(stderr,"%s\n",p);
	MPI_Abort(world,status);
}
