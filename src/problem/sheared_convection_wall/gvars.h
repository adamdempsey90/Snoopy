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

#define		NX				128			/**< X Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NY				128			/**< Y Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NZ				128			/**< Z Dimension in real space. */

//#define		MHD						/**< Uncomment to activate MHD*/

#define		BOUSSINESQ				/**< Uncomment to activate Boussinesq */
#define		VERTSTRAT				/**< Vertical stratification. Otherwise, Boussinesq stratifi-cation is in X */

#define		WITH_ROTATION

//#define		BOUNDARY_C

//#define		WITH_SHEAR				/**< Uncomment to activate mean SHEAR */
//#define		TIME_DEPENDANT_SHEAR	/**< Enable Time dependant shear */

//#define		FORCING					/**< Uncomment to use internal forcing of the velocity field (see forcing in timestep.c) */

#define         OUTPUT_TIMEVAR_FILENAME         "data/timevar"
#define         OUTPUT_VTK_FILENAME             "data/v%04i.vtk"
#define	        OUTPUT_SPECTRUM_FILENAME        "data/spectrum.dat"
#define	OUTPUT_DUMP				"dump.dmp"			/**< Dump files filename. */
#define OUTPUT_DUMP_SAV				"dump_sav.dmp"      /**< Previous (saved) output dump. */
#define OUTPUT_DUMP_WRITE			"dump_write.dmp"	/**< dump currently written. */

#define		FFT_PLANNING	FFTW_MEASURE  /**< can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc). Measure leads to longer initialisation of fft routines */

#endif
