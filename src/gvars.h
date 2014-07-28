/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef _GVARS_
#define _GVARS_

#define		NX			     128	        /**< X Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NY			     128		/**< Y Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NZ			     256		/**< Z Dimension in real space. */
//AJB SHOULD BE TWICE DESIRED RESOLUTION IN THE VERTICAL WHEN USING z-WALLS

//#define		MHD						/**< Uncomment to activate MHD*/

#define		BOUSSINESQ				/**< Uncomment to activate Boussinesq */
#define		VERTSTRAT				/**< Vertical stratification. Otherwise, Boussinesq stratification is in X */

//#define		WITH_SHEAR
#define		WITH_ROTATION

#define         BOUNDARY_C

/* #define         OUTPUT_TIMEVAR_FILENAME         "/projects/b1002/adrian/RB/128om100F1Lh1Lz10nu2p5kap2p5/timevar" */
/* #define         OUTPUT_VTK_FILENAME             "/projects/b1002/adrian/RB/128om100F1Lh1Lz10nu2p5kap2p5/v%04i.vtk" */
/* #define	        OUTPUT_SPECTRUM_FILENAME        "/projects/b1002/adrian/RB/128om100F1Lh1Lz10nu2p5kap2p5/spectrum.dat" */
/* #define	        OUTPUT_PROFILE_FILENAME         "/projects/b1002/adrian/RB/128om100F1Lh1Lz10nu2p5kap2p5/profiles.dat" */
/* #define	OUTPUT_DUMP				"/projects/b1002/adrian/RB/128om100F1Lh1Lz10nu2p5kap2p5/dump.dmp" */
/* #define OUTPUT_DUMP_SAV				"/projects/b1002/adrian/RB/128om100F1Lh1Lz10nu2p5kap2p5/dump_sav.dmp" */
/* #define OUTPUT_DUMP_WRITE                       "/projects/b1002/adrian/RB/128om100F1Lh1Lz10nu2p5kap2p5/dump_write.dmp" */

#define         OUTPUT_TIMEVAR_FILENAME         "/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/timevar"
#define         OUTPUT_VTK_FILENAME             "/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/v%04i.vtk"
#define	        OUTPUT_SPECTRUM_FILENAME        "/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/spectrum.dat"
#define	        OUTPUT_HSPECTRUM_FILENAME        "/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/hspectrum.dat"
#define	        OUTPUT_PROFILE_FILENAME        "/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/profiles.dat" //AJB new

#define	OUTPUT_DUMP				"/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/dump.dmp"			/**< Dump files filename. */
#define OUTPUT_DUMP_SAV				"/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/dump_sav.dmp"      /**< Previous (saved) output dump. */
#define OUTPUT_DUMP_WRITE			"/projects/b1002/horizonal_spectrum_restarts/snoopyRBHC/data/dump_write.dmp"	/**< dump currently written. */

//#define		FORCING					/**< Uncomment to use internal forcing of the velocity field (see forcing in timestep.c) */

#define		FFT_PLANNING	FFTW_MEASURE  /**< can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc). Measure leads to longer initialisation of fft routines */
#define FFTW_MPI_SUPPORT


#endif
