#include <stdlib.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include <math.h>
#define MAX_N_BIN					10000
#define OUTPUT_SPECTRUM_K_BIN		(2.0 * M_PI)
//#ifndef OUTPUT_SPECTRUM_FILENAME //defined in gvars.h //AJB
//#define	OUTPUT_SPECTRUM_FILENAME	"spectrum.dat"
//#endif

/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_spectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + creal( 2.0 * wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}

/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/
/* #ifdef WITH_ELLIPTICAL_VORTEX */
/* void output1Dspectrum(const struct Field fldi, const double ti) { */
/* 	int i; */
/* 	// The zero array is to be used for dummy spectrums */
/* 	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) { w1[i] = 0.0; } */
/* 	write_spectrum(fldi.vx, fldi.vx, ti); */
/* 	write_spectrum(fldi.vy, fldi.vy, ti); */
/* 	write_spectrum(fldi.vz, fldi.vz, ti);	 */
/* 	// Transport spectrums */
/* 	write_spectrum(fldi.vx,fldi.vy, ti); */
/* 	write_spectrum(fldi.vx,fldi.vy, ti); */
/* 	// 3D/shell Transfer spectrums */
/* 	// Kinetic energy transfer */
/* 	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) { */
/* 		w1[i] =  fldi.vx[i]; */
/* 		w2[i] =  fldi.vy[i]; */
/* 		w3[i] =  fldi.vz[i]; */
/* 	} */
/* 	gfft_c2r_t(w1);	gfft_c2r_t(w2); gfft_c2r_t(w3); */
/* 	/\* Compute the convolution for the advection process *\/ */
/* 	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) { */
/* 	  wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL); */
/* 	  wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL); */
/* 	  wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL); */
/* 	  wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL); */
/* 	  wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL); */
/* 	  wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL); */
/* 	} */
/* 	gfft_r2c_t(wr4); gfft_r2c_t(wr5); gfft_r2c_t(wr6);  */
/* 	gfft_r2c_t(wr7); gfft_r2c_t(wr8); gfft_r2c_t(wr9); */
/* 	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) { */
/* 		w1[i] = - I * mask[i] * ( */
/* 					kxt[i] * w4[i] + kyt[i] * w7[i] + kzt[i] * w8[i] ); */
/* 		w2[i] = - I * mask[i] * ( */
/* 					kxt[i] * w7[i] + kyt[i] * w5[i] + kzt[i] * w9[i] ); */
/* 		w3[i] = - I * mask[i] * ( */
/* 					kxt[i] * w8[i] + kyt[i] * w9[i] + kzt[i] * w6[i] ); */
/* 	} */
/* 	write_spectrum(fldi.vx, w1, ti); */
/* 	write_spectrum(fldi.vy, w2, ti); */
/* 	write_spectrum(fldi.vz, w3, ti); */
/*         write_spectrum(fldi.vz, w3, ti); //AJB EMPTY FIELDS TO FILL WITH ANISOTROPIC SPECTRA ETC */
/*         write_spectrum(fldi.vz, w3, ti); */
/* 	write_spectrum(fldi.vz, w3, ti); */
/* 	write_spectrum(fldi.vz, w3, ti); */
/* 	write_spectrum(fldi.vz, w3, ti); */
/* 	write_spectrum(fldi.vz, w3, ti); */
/* #ifdef MHD */
/* 	write_spectrum(fldi.bx, fldi.bx, ti); */
/* 	write_spectrum(fldi.by, fldi.by, ti); */
/* 	write_spectrum(fldi.bz, fldi.bz, ti); */
/* 	write_spectrum(fldi.bx, fldi.by, ti); //bxby */
/* #else */
/* 	write_spectrum(w1, w1, ti); */
/* 	write_spectrum(w1, w1, ti); */
/* 	write_spectrum(w1, w1, ti); */
/* 	write_spectrum(w1, w1, ti); */
/* #endif */

/* 	return; */
/* } */
/* #else //not WITH_ELLIPTICAL_VORTEX */
void output1Dspectrum(const struct Field fldi, const double ti) {
	int i;
	
	DEBUG_START_FUNC;
	
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	// V,B and theta spectrums
	write_spectrum(fldi.vx, fldi.vx, ti);
	write_spectrum(fldi.vy, fldi.vy, ti);
	write_spectrum(fldi.vz, fldi.vz, ti);
	
#ifdef MHD
	write_spectrum(fldi.bx, fldi.bx, ti);
	write_spectrum(fldi.by, fldi.by, ti);
	write_spectrum(fldi.bz, fldi.bz, ti);
#else
#ifdef BOUSSINESQ //AJB 10/04/13
        write_spectrum(fldi.vx, fldi.th, ti);
        write_spectrum(fldi.vy, fldi.th, ti);
        write_spectrum(fldi.vz, fldi.th, ti);
#else  //AJB 10/04/13
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif  //AJB 10/04/13
#endif

#ifdef BOUSSINESQ
	write_spectrum(fldi.th, fldi.th, ti);
#else 
	write_spectrum(w1, w1, ti);
#endif

	// Transport spectrums
	write_spectrum(fldi.vx,fldi.vy, ti);
#ifdef MHD
	write_spectrum(fldi.bx,fldi.by, ti);
#else
	write_spectrum(w1, w1, ti);
#endif
	// Transfer spectrums
	// Kinetic energy transfer
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w4[i] + kyt[i] * w7[i] + kzt[i] * w8[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w7[i] + kyt[i] * w5[i] + kzt[i] * w9[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w8[i] + kyt[i] * w9[i] + kzt[i] * w6[i] );
	}
	
	write_spectrum(fldi.vx, w1, ti);
	write_spectrum(fldi.vy, w2, ti);
	write_spectrum(fldi.vz, w3, ti);
	
#ifdef MHD

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// (vx,vy,vz) is in w1-w3 and (bx,by,bz) is in (w4-w6). It is now time to compute the emfs in w7-w9...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = (wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL*NTOTAL);
		wr8[i] = (wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL*NTOTAL);
		wr9[i] = (wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf involved in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * mask[i] * (kyt[i] * w9[i] - kzt[i] * w8[i]);
		w2[i] = I * mask[i] * (kzt[i] * w7[i] - kxt[i]* w9[i]);
		w3[i] = I * mask[i] * (kxt[i]* w8[i] - kyt[i] * w7[i]);
	}

	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	
	// Let's do the Lorentz Force
	// We already have (bx,by,bz) in w4-w6. No need to compute them again...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr4[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr2[i] = wr5[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr3[i] = wr6[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr4[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr4[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr5[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}

	gfft_r2c_t(wr1);
	gfft_r2c_t(wr2);
	gfft_r2c_t(wr3);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);


	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = I * mask[i] * (kxt[i] * w1[i] + kyt[i] * w7[i] + kzt[i] * w8[i]);
		w5[i] = I * mask[i] * (kxt[i] * w7[i] + kyt[i] * w2[i] + kzt[i] * w9[i]);
		w6[i] = I * mask[i] * (kxt[i] * w8[i] + kyt[i] * w9[i] + kzt[i] * w3[i]);
	}
	
	write_spectrum(fldi.vx, w4, ti);
	write_spectrum(fldi.vy, w5, ti);
	write_spectrum(fldi.vz, w6, ti);
	
	// Helicity spectrums
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * ik2t[i] * (kyt[i] * fldi.bz[i] - kzt[i] * fldi.by[i] );
		w2[i] = I * ik2t[i] * (kzt[i] * fldi.bx[i] - kxt[i] * fldi.bz[i] );
		w3[i] = I * ik2t[i] * (kxt[i] * fldi.by[i] - kyt[i] * fldi.bx[i] );
	}
	
	write_spectrum(fldi.bx, w1, ti);
	write_spectrum(fldi.by, w2, ti);
	write_spectrum(fldi.bz, w3, ti);
	//	write_uAu_spectrum(fldi.bx, fldi.by, ti); //BAB
#else //MHD
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);

#endif //MHD
	
	DEBUG_END_FUNC;
	
	return;
}
/* #endif // WITH_ELLIPTICAL_VORTEX */
/**********************************************************/
/**
	Initialise the 1D spectrum output routine, used
	to output the spectrum
	This routine print the mode ks in the first line
	It also counts the number of mode in each shell and 
	output it in the second line of OUTPUT_SPECTRUM_FILENAME
*/
/*********************************************************/
void init1Dspectrum() {
	int i,j,k,m;
	int nbin;
	FILE * ht;
	double spectrum[ MAX_N_BIN ];
	
	DEBUG_START_FUNC;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < nbin; m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");
	}
	
	for( i = 0; i < MAX_N_BIN ; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC ; i++) {
	  for( j = 0; j < NY_COMPLEX; j++) {
	    for( k = 0; k < NZ_COMPLEX; k++) {
	      m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
	      if ( m < nbin)
		spectrum[ m ] = spectrum[ m ] + 1.0;
	    }
	  }
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin ; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		for( i = 0; i < nbin ; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
	
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	
	DEBUG_END_FUNC;
	
	return;
}

//AJB 12/06/13
/***********************************************************/
/**
	Compute 1D horizontally averaged profiles of uz,th etc and output to OUTPUT_PROFILE_FILENAME
	@param wi Input field to average and output
	@param ti Current time
*/
/***********************************************************/
void init1Dprofile() {
	int m;
	FILE * ht;
	if(rank==0) {
	  ht = fopen(OUTPUT_PROFILE_FILENAME,"w");
	}
}
void output1Dprofile(const struct Field fldi, const double ti) {
  int i,j,k,indx;
  double z;
  double pr[28][NZ];
  FILE *ht;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(static)
#endif
  for(i=0;i<NTOTAL_COMPLEX;i++){
    w1[i]=fldi.vx[i]; //ux
    w2[i]=fldi.vy[i]; //uy
    w3[i]=fldi.vz[i]; //uz
    w4[i]=fldi.th[i]; //th
    w5[i]=I*kzt[i]*fldi.th[i]; //dzth
    w6[i]=I*kxt[i]*fldi.vx[i]; //dxux
    w7[i]=I*kyt[i]*fldi.vx[i]; //dyux
    w8[i]=I*kzt[i]*fldi.vx[i]; //dzux
    w9[i]=I*kxt[i]*fldi.vy[i]; //dxuy
    w10[i]=I*kyt[i]*fldi.vy[i]; //dyuy
    w11[i]=I*kzt[i]*fldi.vy[i]; //dzuy
    w12[i]=I*kxt[i]*fldi.vz[i]; //dxuz
    w13[i]=I*kyt[i]*fldi.vz[i]; //dyuz
    w14[i]=I*kzt[i]*fldi.vz[i]; //dzuz
    //don't use w15 as heating/cooling uses it.
    w16[i]=I*kxt[i]*fldi.th[i]; //dxth
    w17[i]=I*kyt[i]*fldi.th[i]; //dyth
  }
  gfft_c2r_t(w1);gfft_c2r_t(w2);gfft_c2r_t(w3);gfft_c2r_t(w4);gfft_c2r_t(w5);
  gfft_c2r_t(w6);gfft_c2r_t(w7);gfft_c2r_t(w8);gfft_c2r_t(w9);gfft_c2r_t(w10);
  gfft_c2r_t(w11);gfft_c2r_t(w12);gfft_c2r_t(w13);gfft_c2r_t(w14);gfft_c2r_t(w16);gfft_c2r_t(w17);
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(static)
#endif
  for(i=0;i<2*NTOTAL_COMPLEX;i++){
    wr1[i]/=((double) NTOTAL ); //ux
    wr2[i]/=((double) NTOTAL ); //uy
    wr3[i]/=((double) NTOTAL ); //uz
    wr4[i]/=((double) NTOTAL ); //th
    wr5[i]/=((double) NTOTAL ); //dzth
    wr6[i]/=((double) NTOTAL ); //dxux
    wr7[i]/=((double) NTOTAL ); //dyux
    wr8[i]/=((double) NTOTAL ); //dzux
    wr9[i]/=((double) NTOTAL ); //dxuy
    wr10[i]/=((double) NTOTAL ); //dyuy
    wr11[i]/=((double) NTOTAL ); //dzuy
    wr12[i]/=((double) NTOTAL ); //dxuz
    wr13[i]/=((double) NTOTAL ); //dyuz
    wr14[i]/=((double) NTOTAL ); //dzuz
    wr16[i]/=((double) NTOTAL ); //dxth
    wr17[i]/=((double) NTOTAL ); //dyth
  }
  for(i=0;i<28;i++){
      for(k=0;k<NZ;k++){
	pr[i][k]=0.0;}}
  for(i=0;i<NX/NPROC;i++){//Horizontal averages
    for(j=0;j<NY;j++){
      for(k=0;k<NZ;k++){
	z=(param.lz*k)/NZ;
	indx=k+(NZ+2)*j+(NZ+2)*NY*i;
	pr[0][k]+=wr1[indx]*wr1[indx]/((double) NX*NY); //<ux^2>
	pr[1][k]+=wr2[indx]*wr2[indx]/((double) NX*NY); //<uy^2>
	pr[2][k]+=wr3[indx]*wr3[indx]/((double) NX*NY); //<uz^2>
	pr[3][k]+=wr4[indx]*wr4[indx]/((double) NX*NY); //<th^2>
	pr[4][k]+=(param.N2*z+wr4[indx])/((double) NX*NY); //<thtot>
	pr[5][k]+=wr4[indx]*wr3[indx]/((double) NX*NY); //<th uz>
	pr[6][k]+=wr4[indx]/((double) NX*NY); //<th>
	pr[7][k]+=wr1[indx]/((double) NX*NY); //<ux>
	pr[8][k]+=wr2[indx]/((double) NX*NY); //<uy>
	pr[9][k]+=wr5[indx]/((double) NX*NY); //<dzth>
	pr[10][k]+=wr6[indx]*wr6[indx]/((double) NX*NY); //<dxux^2>
	pr[11][k]+=wr10[indx]*wr10[indx]/((double) NX*NY); //<dyuy^2>
	pr[12][k]+=wr14[indx]*wr14[indx]/((double) NX*NY); //<dzuz^2>
	pr[13][k]+=wr7[indx]*wr7[indx]/((double) NX*NY); //<dyux^2>
	pr[14][k]+=wr8[indx]*wr8[indx]/((double) NX*NY); //<dzux^2>
	pr[15][k]+=wr9[indx]*wr9[indx]/((double) NX*NY); //<dxuy^2>
	pr[16][k]+=wr11[indx]*wr11[indx]/((double) NX*NY); //<dzuy^2>
	pr[17][k]+=wr12[indx]*wr12[indx]/((double) NX*NY); //<dxuz^2>
	pr[18][k]+=wr13[indx]*wr13[indx]/((double) NX*NY); //<dyuz^2>
	pr[19][k]+=wr16[indx]*wr16[indx]/((double) NX*NY); //<dxth^2>
	pr[20][k]+=wr17[indx]*wr17[indx]/((double) NX*NY); //<dyth^2>
	pr[21][k]+=wr5[indx]*wr5[indx]/((double) NX*NY); //<dzth^2>
	pr[22][k]+=wr7[indx]*wr9[indx]/((double) NX*NY); //<dxuydyux>
	pr[23][k]+=wr8[indx]*wr12[indx]/((double) NX*NY); //<dxuzdzux>
	pr[24][k]+=wr13[indx]*wr11[indx]/((double) NX*NY); //<dyuzdzuy>
	pr[25][k]+=wr1[indx]*wr2[indx]/((double) NX*NY); //<uxuy>
	pr[26][k]+=wr1[indx]*wr3[indx]/((double) NX*NY); //<uxuz>
	pr[27][k]+=wr2[indx]*wr3[indx]/((double) NX*NY); //<uyuz>
      }
    }
  }
#ifdef MPI_SUPPORT
	for(i=0;i<28;i++){
	  for(k=0;k<NZ;k++) {//could be made more efficient...
	    reduce(&pr[i][k],1);}}
#endif
  //Now begin outputting 1D profiles to file
  if(rank==0) {
    ht = fopen(OUTPUT_PROFILE_FILENAME,"a");
    for(i=0;i<NZ;i++) {
      z=(param.lz*i)/NZ;
      fprintf(ht,"%08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t %08e\t \n",z,pr[0][i],pr[1][i],pr[2][i],pr[3][i],pr[4][i],pr[5][i],pr[6][i],pr[7][i],pr[8][i],pr[9][i]+param.N2,pr[10][i],pr[11][i],pr[12][i],pr[13][i],pr[14][i],pr[15][i],pr[16][i],pr[17][i],pr[18][i],pr[19][i],pr[20][i],pr[21][i],pr[22][i],pr[23][i],pr[24][i],pr[25][i],pr[26][i],pr[27][i]);
    }
    fclose(ht);
  }
  return;
} //AJB

void initHSpectrum() {
	int m;
	FILE * hs;
	if(rank==0) {
	  hs = fopen(OUTPUT_HSPECTRUM_FILENAME,"w");
	}
	return;
}
void write_hspectrum(const double complex wi[],const double complex wj[], const double ti, const double lowz, const double hiz) {
	int i,j,k,m,indx,indx2, indx1, zcount, nxny, final_count, nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *hs;
	double khmax,znow;
	double complex *tempc1, *tempc2;
	double *tempr1, *tempr2, *kx2, *ky2, *kh2, *final_whr1;
	double *temp1, *temp2;
	double local_kx,local_ky,local_kh;	
	
// Allocate arrays.	
	nxny= NX*(NY/2+1);
	tempc1 = (double complex *) fftw_malloc( sizeof(double complex) * nxny);
	tempr1 = (double *) tempc1;
	tempc2 = (double complex *) fftw_malloc( sizeof(double complex) * nxny);
	tempr2 = (double *) tempc2;
	kx2 = (double *)fftw_malloc(sizeof(double)*nxny);
	ky2 = (double *)fftw_malloc(sizeof(double)*nxny);
	kh2 = (double *)fftw_malloc(sizeof(double)*nxny);
	final_whr1 = (double *)fftw_malloc(sizeof(double)*nxny);
	final_count=0;

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for(i=0;i<2*nxny; i++) {	// Initialize
		tempr1[i]=0; 
		tempr2[i]=0;
		if(i<nxny) {
			wrh3[i]=0; 
			final_whr1[i]=0;
		}
	}
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif 
	for(i=0;i<NX*2*(NY/2+1)*NZ;i++) {
		wrh4[i]=0;
		wrh5[i]=0;
	}

// Copy input constant arrays into work arrays wh1 and wh2
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for(i=0;i<NTOTAL_COMPLEX;i++) {
		wh1[i] = wi[i];
		wh2[i] = wj[i];
	}
	
	

// First transform arrays to real space, no transposing here.

	gfft_c2r(wh1); gfft_c2r(wh2);
		

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for(i=0;i<2*NTOTAL_COMPLEX; i++) {
 		wrh1[i]/=((double) NTOTAL);
 		wrh2[i]/= ((double) NTOTAL);
	}
// Now have arrays that are x by y by z, distributed along x.
// Transpose the arrays so that the z-dimension is distributed 
// And x-y matrices are contiguous in memory. So arrays are now z by x by y.

#ifdef MPI_SUPPORT
	transpose_real_XZ(wrh1,wrh4); transpose_real_XZ(wrh2,wrh5);		// use mpi transposes
#else													// else just do it.
	temp1=(double *)fftw_malloc(sizeof(double)*2*NX*(NY/2+1)*NZ);
	temp2=(double *)fftw_malloc(sizeof(double)*2*NX*(NY/2+1)*NZ);
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)
#endif
	for(i=0;i<2*NX*(NY/2+1)*NZ;i++) {
		temp1[i]=0;
		temp2[i]=0;
	}
	
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)
#endif
	for(i=0;i<NX_COMPLEX;i++) {
		for(j=0;j<NY_COMPLEX;j++) {
			for(k=0;k<NZ;k++) {
				temp1[j+2*(NY/2+1)*i + (NY/2+1)*2*NZ*k] = wrh1[i*NY_COMPLEX*2*NZ_COMPLEX+j*NZ_COMPLEX+k];
				temp2[j+2*(NY/2+1)*i + (NY/2+1)*2*NZ*k] = wrh2[i*NY_COMPLEX*2*NZ_COMPLEX+j*NZ_COMPLEX+k];
			}
		}
	}
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for(i=0;i<2*NX*(NY/2+1)*NZ;i++) {
		wrh4[i] = temp1[i];
		wrh5[i] = temp2[i];
	}
	fftw_free(temp1); fftw_free(temp2);
#endif

	zcount=0;
	
// Loop over the local z-indicies, compute the height that we're at and then if it is 
// within the bounds given then do the 2D transform in the horizontal directions

	for(k=0; k < NZ/NPROC; k++) {
	
#ifdef BOUNDARY_C
		znow = param.lz *( k + rank*NZ/NPROC)/NZ;
#else
		 znow=- param.lz / 2 + (param.lz * (k + rank * NZ / NPROC)) / NZ;
			
#endif
		if(znow>= lowz && znow <= hiz) {	
			zcount++;
			for(i=0; i<NX; i++) {		// First make wri*wrj array with padding in y direction.
				for(j=0; j<2*(NY/2+1); j++) {
					indx1=j + 2*(NY/2+1)*i + NX*2*(NY/2+1)*k;
					tempr1[j+2*(NY/2+1)*i] = wrh4[indx1];
					tempr2[ j+2*(NY/2+1)*i] = wrh5[indx1]; 
				}
			}

			gfft_r2c_2Dslice(tempr1);		// Transform
			gfft_r2c_2Dslice(tempr2);
			
			
			for(i=0;i<nxny;i++) {	// Copy over |FT(wri*wrj)|^2
				wrh3[i]+=(double)creal(tempc1[i]*conj(tempc2[i])); 
					
			}
			


			
		}
	}
// Send to root and average over number of planes in the bounds.
#ifdef MPI_SUPPORT		
	MPI_Reduce(wrh3,final_whr1,NX*(NY/2+1),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&zcount,&final_count,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
#endif
	
	if(rank==0) {
// Output final spectrum
		khmax=2*M_PI*pow((NX/2-1)*(NX/2-1)/(param.lx*param.lx)+(NY/2-1)*(NY/2-1)/(param.ly*param.ly),.5);
#ifdef _OPENMP
		#pragma omp parallel for private(i) schedule(static)
#endif
		for(i=0;i<nxny;i++) {
			final_whr1[i] /= ((double)final_count);
		}
		
		for(i=0;i<NX; i++) {
			for(j=0; j<(NY/2+1); j++) {
				kx2[j+(NY/2+1)*i] =(2.0 * M_PI) / param.lx*(fmod( i + (NX / 2) ,  NX ) - NX / 2 );
				ky2[j+(NY/2+1)*i] =(2.0 * M_PI) / param.ly * j;
				kh2[j+(NY/2+1)*i] = pow(kx2[j+(NY/2+1)*i]*kx2[j+(NY/2+1)*i]+ky2[j+(NY/2+1)*i]*ky2[j+(NY/2+1)*i],.5);
			}
		}
			
//		nbin = (int) ceil(khmax / (2*M_PI));
		nbin = (int) ceil(pow((NX/2-1)*(NX/2-1)+(NY/2-1)*(NY/2-1),.5));	
		for( i = 0; i < 10000; i++ )	spectrum[ i ] = 0.0;
	
		for(i=0;i<NX; i++) {
			for(j=0;j<(NY/2+1);j++) {
				indx1=j+i*(NY/2+1);
				local_kx = (double)(fmod(i+NX/2,NX)-NX/2);
				local_ky = (double)j;
				local_kh = pow(local_kx*local_kx+local_ky*local_ky,.5);
				m = (int) floor(local_kh + 0.5 );
		
				if( j == 0) {
					// j=0, we have all the modes.
					spectrum[ m ] +=  final_whr1[indx1]/((double)NX*NY*NX*NY);
				}
				else {
								// j>0, only half of the complex plane is represented.
					spectrum[ m ] += 2.0*final_whr1[indx1]/((double)NX*NY*NX*NY);
				}
			}
		}
		
		hs = fopen(OUTPUT_HSPECTRUM_FILENAME,"a");
		fprintf(hs,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(hs,"%08e\t",spectrum[i]);
	
		fprintf(hs,"\n");
		
		if(ferror(hs)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing horizontal spectrum file");
		
		fclose(hs);
	}	

	fftw_free(kx2); fftw_free(ky2); fftw_free(kh2); fftw_free(final_whr1); 
	fftw_free(tempr1); fftw_free(tempr2);
	
	return;
}



void outputHSpectrum(const struct Field fldi, const double ti) {
 	int i;
 	FILE *hs;

	MPI_Printf("Outputting Horizontal Spectrum at time %g\n", ti);

//param.delz,(param.lz/2.0)-param.delz
 
// 	write_HSpectrum(fldi.vx,fldi.vx,ti,-.25,.25);
// 	write_HSpectrum(fldi.vy,fldi.vy,ti,-.25,.25);
// 	write_HSpectrum(fldi.vz,fldi.vz,ti,-.25,.25);
// 	write_HSpectrum(fldi.vz,fldi.th,ti,-.25,.25);
// 	write_HSpectrum(w5,w5,ti,-.25,.25);
	
 	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	// V,B and theta spectrums
	write_hspectrum(fldi.vx, fldi.vx, ti,param.delz,(param.lz/2.0)-param.delz);
	write_hspectrum(fldi.vy, fldi.vy, ti,param.delz,(param.lz/2.0)-param.delz);
	write_hspectrum(fldi.vz, fldi.vz, ti,param.delz,(param.lz/2.0)-param.delz);
	

#ifdef BOUSSINESQ //AJB 10/04/13
        write_hspectrum(fldi.vx, fldi.th, ti,param.delz,(param.lz/2.0)-param.delz);
        write_hspectrum(fldi.vy, fldi.th, ti,param.delz,(param.lz/2.0)-param.delz);
        write_hspectrum(fldi.vz, fldi.th, ti,param.delz,(param.lz/2.0)-param.delz);
#else  //AJB 10/04/13
	write_hspectrum(w1, w1, ti);
	write_hspectrum(w1, w1, ti);
	write_hspectrum(w1, w1, ti);
#endif  //AJB 10/04/13

#ifdef BOUSSINESQ
	write_hspectrum(fldi.th, fldi.th, ti,param.delz,(param.lz/2.0)-param.delz);
#else 
	write_hspectrum(w1, w1, ti);
#endif

	// Transport spectrums
	write_hspectrum(fldi.vx,fldi.vy, ti,param.delz,(param.lz/2.0)-param.delz);

	// Transfer spectrums
	// Kinetic energy transfer
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w4[i] + kyt[i] * w7[i] + kzt[i] * w8[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w7[i] + kyt[i] * w5[i] + kzt[i] * w9[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w8[i] + kyt[i] * w9[i] + kzt[i] * w6[i] );
	}
	
	write_hspectrum(fldi.vx, w1, ti,param.delz,(param.lz/2.0)-param.delz);
	write_hspectrum(fldi.vy, w2, ti,param.delz,(param.lz/2.0)-param.delz);
	write_hspectrum(fldi.vz, w3, ti,param.delz,(param.lz/2.0)-param.delz);
	


	return;
}
