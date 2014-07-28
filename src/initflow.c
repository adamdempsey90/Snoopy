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
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "common.h"
#include "gfft.h"
#include "output/output_dump.h"
#include "symmetries.h"

#include "debug.h"

/** Allow one to init a structure in real space using ordinary defined x,y,z coordinates */

void init_SpatialStructure(struct Field fldi) {
        double *x,*y,*z;
	int i,j,k;
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Allocate coordinate arrays
	x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");

	// Initialize the arrays
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.ly / 2 + (param.ly * j ) / NY;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lz / 2 + (param.lz * k ) / NZ;
#ifdef BOUNDARY_C
				z[k+(NZ+2)*j+(NZ+2)*NY*i] = param.lz*k/NZ;
//				if(k<NZ/2) {//z=0...1 AJB 18/06/13
//				  z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = (param.lz/2.0)*k/NZ/2.0; } else {
//				  z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
//			}
#endif
			}
		}
	}
	
	// Initialize the extra points (k=NZ and k=NZ+1) to zero to prevent stupid things from happening...
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = NZ ; k < NZ + 2 ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	// Init work array to zero
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ + 2 ; k++) {
				wr1[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr2[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr3[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr4[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr5[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr6[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr7[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	/*******************************************************************
	** This part can be modified              **************************
	********************************************************************/
	
	// The velocity field vx,vy,vz is stored in wr1,wr2,wr3
	// The magnetic field bx,by,bz is stored in wr4,wr5,wr6 (ignored if MHD is not set)

//	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		// Example: init a flux tube in the x direction+a vertical displacement
	  //	wr4[i] = exp(-(y[i]*y[i]+z[i]*z[i])*20.0);
	  //	wr3[i] = 0.5*cos(x[i] * 2.0 * M_PI);
//			wr3[i] = sin(2*M_PI*x[i]/param.lx)*sin(2*M_PI*y[i]/param.ly)*sin(2*M_PI*z[i]/param.lz);
// //AJB 17/06/13 TEMP PROFILE PIECEWISE LINEAR
// 	   if(z[i]<0.25){//z<0.25 
// 	     wr7[i]=-1.0*z[i]; 
// 	   }else if(z[i]>=0.25 && z[i]<0.75){//z>0.25->0.75 
// 	    wr7[i]=-0.25+(z[i]-0.25); 
// 	   }else{//z>0.75 
// 	     wr7[i]=0.25-1.0*(z[i]-0.75); 
// 	  } 
// 	    wr7[i]*=2.0; //account for defn of sin 

	  /* //Taylor-Green vortex appropriate for periodic BCs. */
	  /* wr1[i] = -psi0*by*cos(bx*(x[i]-dx))*sin(by*(y[i]-dy)); */
	  /* wr2[i] = psi0*bx*sin(bx*(x[i]-dx))*cos(by*(y[i]-dy)); */

		// Example: twisted flux tube + vertical displacement
		/* wr4[i] = exp(-(y[i]*y[i]+z[i]*z[i])/(0.2*0.2)); */
		/* wr5[i] = fabs(z[i])*1.0*wr4[i]; */
		/* wr6[i] = -fabs(y[i])*1.0*wr4[i]; */
		/* wr3[i] = 0.5*cos(x[i] * 2.0 * M_PI); */
		/* if (i==3*NY*(NZ+2)) fprintf(stderr," %d %e",i,x[i]); */
//	}

	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Fourier transform everything
	gfft_r2c(wr1);
	gfft_r2c(wr2);
	gfft_r2c(wr3);
	gfft_r2c(wr4);
	gfft_r2c(wr5);
	gfft_r2c(wr6);
	gfft_r2c(wr7);

	// Transfer data in the relevant array (including dealiasing mask)
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
	  fldi.vx[i] += w1[i] * mask[i];
	  fldi.vy[i] += w2[i] * mask[i];
	  fldi.vz[i] += w3[i] * mask[i];
#ifdef BOUSSINESQ
	  fldi.th[i] += w7[i] * mask[i];
#endif
#ifdef MHD
		fldi.bx[i] += w4[i] * mask[i];
		fldi.by[i] += w5[i] * mask[i];
		fldi.bz[i] += w6[i] * mask[i];
#endif
	}

	// free memory
	fftw_free(x);
	fftw_free(y);
	fftw_free(z);
	
	//done
	return;
}

#ifdef WITH_ELLIPTICAL_VORTEX //AJB 
void init_Kelvin_Wave(struct Field fldi) {
  double *x,*y,*z,theta,Phi,k0,factor,amp;
  int i,j,k,indx;
  x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
  if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
  y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
  if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
  z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
  if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");
  for(i = 0 ; i < NX/NPROC ; i++) {
    for(j = 0 ; j < NY ; j++) {
      for(k = 0 ; k < NZ ; k++) {
	x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
	y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.ly / 2 + (param.ly * j ) / NY;
	z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lz / 2 + (param.lz * k ) / NZ;
      }
    }
  }
  for(i = 0 ; i < NX/NPROC ; i++) {
    for(j = 0 ; j < NY ; j++) {
      for(k = NZ ; k < NZ + 2 ; k++) {
	x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
	y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
	z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
      }
    }
  }
  // INPUT Kelvin Wave in Fourier space.
  //IDX3D = (k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i);
  //theta = M_PI/4.0; Phi = 0.0*M_PI/4.0; k0 = 4.0*(2.0*M_PI/param.lx); 
  i = param.Kelvin_Wave_i; j = param.Kelvin_Wave_j; k = param.Kelvin_Wave_k;
  /* MPI_Printf("ijk: %i \t %i \t %i\n",i,j,k); */
  theta = atan((param.lz/param.lx)*(( (double) i)/((double) k)));
  k0 = 2.0*M_PI*sqrt((((double)(i*i))/(param.lx*param.lx))+(((double)(j*j))/(param.ly*param.ly))+(((double)(k*k))/(param.lz*param.lz)));
  Phi = atan((param.lx/param.ly)*(( (double) j)/((double) i)));
/*   MPI_Printf("theta: %f \n",theta); */
/*   MPI_Printf("Phi: %f \n",Phi); */
/*   MPI_Printf("k0: %f \n",k0); */

  amp = param.Kelvin_Wave_amp;
  factor = 0.5*NTOTAL;
  indx = (k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i);

  fldi.vx[indx] = amp*(((double) k)*(2.0*M_PI/param.lz)/k0)*factor*cos(Phi); //costh factor
  fldi.vy[indx] = amp*sin(Phi)*factor; //note minus sign here... (chosen artificially to make -Pi/4 chosen phase)
  fldi.vz[indx] = -(((double) i)*amp*(2.0*M_PI/param.lx)/k0)*factor*cos(Phi); //sinth factor

  projector(fldi.vx,fldi.vy,fldi.vz);

  return;

}
#endif

void init_B_field(struct Field fldi) { //AJB
#ifdef MHD
  int i,j,k,indx;

  if(rank==0) {
    i = 1; j = 0; k = 0; //single Fourier mode B field
    indx = (k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i);
    
    fldi.bx[indx] = 0.0;
    fldi.by[indx] = 0.0;
    fldi.bz[indx] = param.B_amp*((double) NTOTAL);
  }
  
  projector(fldi.bx,fldi.by,fldi.bz);
#endif
}

void init_Lamb_Vortex(struct Field fldi) { //AJB
  double x,y,z,a,b;
  int i,j,k;

  a = param.Lamb_Vortex_aspect*param.Lamb_Vortex_lambda; b = param.Lamb_Vortex_lambda;
  for(i = 0 ; i < NX/NPROC ; i++) {
    x = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
    for(j = 0 ; j < NY ; j++) {
      y = - param.ly / 2 + (param.ly * j) / NY;
      for(k = 0 ; k < NZ ; k++) {
	//z = - param.lz / 2 + (param.lz * k ) / NZ;
	if(x * x / (a * a) + y * y / (b * b) < 1) {
	  wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = param.Lamb_Vortex_vorticity;
	}
	else {
	  wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = 0.0;
	}
      }
    }
  }
  // transform
  gfft_r2c_t(wr1);
  for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
    fldi.vx[ i ] +=  I * ky[i] * w1[i] * ik2t[i];
    fldi.vy[ i ] += -I * kxt[i] * w1[i] * ik2t[i];
  }

  return;
}

void init_KidaVortex(struct Field fldi) {
	double a = param.vortex_a;
	double b = param.vortex_b;
	
	int i,j,k;
	
	double w0, x, y;
	double chi;
	
	chi = b / a;
	w0 = 1.0/chi*(chi + 1.0)/(chi-1.0);			// According to Kida!
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		x = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
		for(j = 0 ; j < NY ; j++) {
			y = - param.ly / 2 + (param.ly * j) / NY;
#ifdef WITH_2D
			if(x * x / (a * a) + y * y / (b * b) < 1) {
					// we are in the vortex
					wr1[j + (NY+2) * i] = -w0;
			}
			else {
				wr1[j + (NY+2) * i] = 0.0;
			}
#else
			for(k = 0 ; k < NZ ; k++) {
				if(x * x / (a * a) + y * y / (b * b) < 1) {
					// we are in the vortex
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = -w0;
				}
				else {
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = 0.0;
				}
			}
#endif
		}
	}
	
	// transform
	gfft_r2c_t(wr1);
	
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[ i ] +=  I * ky[i] * w1[i] * ik2t[i];
		fldi.vy[ i ] += -I * kxt[i] * w1[i] * ik2t[i];
	}
	
	// done
	return;
}

/************************************/
/** Init some crazy structure involving
/** A kida vortex and a vertical structure
/** for the field */
/***********************************/
void init_Bench(struct Field fldi) {
	const double a = 0.3;
	const double b = 0.4;
	
	int i,j,k;
	
	double w0, x, y;
	double chi;
	
	chi = b / a;
	w0 = 1.0/chi*(chi + 1)/(chi-1.0);			// According to Kida!
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		x = - param.lx / 2. + (param.lx * (i + rank * NX / NPROC)) / NX;
		for(j = 0 ; j < NY ; j++) {
			y = - param.ly / 2. + (param.ly * j) / NY;
			for(k = 0 ; k < NZ ; k++) {
				if(x * x / (a * a) + y * y / (b * b) < 1) {
					// we are in the vortex
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = -w0;
				}
				else {
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = 0.0;
				}
			}
		}
	}
	
	// transform
	gfft_r2c(wr1);
	
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[ i ] +=  I * ky[i] * w1[i] * ik2t[i];
		fldi.vy[ i ] += -I * kxt[i] * w1[i] * ik2t[i];
	}
	
	// Brake vertical symmetry
	if(rank==0) {
		fldi.vx[1] = 1000.0 / NTOTAL;
		fldi.vy[1] = 1000.0 / NTOTAL;
#ifdef MHD
		fldi.bx[1] = 1000.0 / NTOTAL;
		fldi.by[1] = 1000.0 / NTOTAL;
#endif
	}
	// done
	return;
}


void init_LargeScaleNoise(struct Field fldi) {
	int i,j,k;
	int num_force=0;
	int total_num_force;
	double fact;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				if(pow(k2t[ IDX3D ], 0.5) / ( 2.0*M_PI ) < 1.0 / param.noise_cut_length && pow(k2t[ IDX3D ], 0.5) / ( 2.0*M_PI ) > 1.0 / param.noise_cut_length_max && k2t[ IDX3D ] != 0 ) {
					fldi.vx[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vy[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vz[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#ifdef MHD
					fldi.bx[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.by[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.bz[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#endif
					if(mask[IDX3D] > 0) num_force++;
				}
			}
		}
	}
	
	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	
	// Divide by the total number of modes
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[ IDX3D ] = fldi.vx[ IDX3D ] / fact;
				fldi.vy[ IDX3D ] = fldi.vy[ IDX3D ] / fact;
				fldi.vz[ IDX3D ] = fldi.vz[ IDX3D ] / fact;
#ifdef MHD
				fldi.bx[ IDX3D ] = fldi.bx[ IDX3D ] / fact;
				fldi.by[ IDX3D ] = fldi.by[ IDX3D ] / fact;
				fldi.bz[ IDX3D ] = fldi.bz[ IDX3D ] / fact;
#endif
			}
		}
	}
	
  enforce_complex_symm(fldi);  
}

void init_LargeScaleNoise_k(struct Field fldi) {
  int i,j,k,indx;
	int num_force=0;
	int total_num_force;
	double fact;//x,y,z;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
			  if( (pow(k2t[ IDX3D ], 0.5)*(param.lx/(2.0*M_PI)) <= param.kmax) && (pow(k2t[ IDX3D ], 0.5)*(param.lx/(2.0*M_PI)) >= param.kmin) && (k2t[ IDX3D ] != 0) ) {
					fldi.vx[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vy[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vz[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
/* #ifdef MHD */
/* 					fldi.bx[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL; */
/* 					fldi.by[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL; */
/* 					fldi.bz[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL; */
/* #endif */
					if(mask[IDX3D] > 0) num_force++;
				}
			}
		}
	}
	
	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	
	// Divide by the total number of modes
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[ IDX3D ] = fldi.vx[ IDX3D ] / fact;
				fldi.vy[ IDX3D ] = fldi.vy[ IDX3D ] / fact;
				fldi.vz[ IDX3D ] = fldi.vz[ IDX3D ] / fact;
/* #ifdef MHD */
/* 				fldi.bx[ IDX3D ] = fldi.bx[ IDX3D ] / fact; */
/* 				fldi.by[ IDX3D ] = fldi.by[ IDX3D ] / fact; */
/* 				fldi.bz[ IDX3D ] = fldi.bz[ IDX3D ] / fact; */
/* #endif */
			}
		}
	}
	
  enforce_complex_symm(fldi);
}

/******************************************
** Large scale 2D (x,y) noise *************
*******************************************/

void init_LargeScale2DNoise(struct Field fldi) {
	int i,j,k;
	int num_force=0;
	int total_num_force;
	double fact;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			k=0;
			if(kz[ IDX3D ] == 0.0) {
				if(pow(k2t[ IDX3D ], 0.5) / ( 2.0*M_PI ) < 1.0 / param.noise_cut_length_2D) {
					fldi.vx[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vy[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#ifdef MHD
					fldi.bx[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.by[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#endif
					if(mask[IDX3D] > 0) num_force++;
				}
			}
		}
	}
	
	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	
	// Divide by the total number of modes
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			k=0;
			if(kz[ IDX3D ] == 0.0) {
				fldi.vx[ IDX3D ] = fldi.vx[ IDX3D ] / fact;
				fldi.vy[ IDX3D ] = fldi.vy[ IDX3D ] / fact;
#ifdef MHD
				fldi.bx[ IDX3D ] = fldi.bx[ IDX3D ] / fact;
				fldi.by[ IDX3D ] = fldi.by[ IDX3D ] / fact;
#endif
			}
		}
	}
	
  enforce_complex_symm(fldi);  
}


void init_WhiteNoise(struct Field fldi) {
	int i,j,k;
	double fact;
	
	// Excite (2/3)^3*NTOTAL modes
	fact = pow(27.0/8.0*NTOTAL, 0.5);
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fldi.vy[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fldi.vz[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
#ifdef MHD
				fldi.bx[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fldi.by[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
				fldi.bz[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
#endif
			}
		}
	}
	
	enforce_complex_symm(fldi);
  
}

void init_MeanField(struct Field fldi) {
#ifdef MHD
	if(rank==0) {
		fldi.bx[0] = param.bx0 * ((double) NTOTAL);
		fldi.by[0] = param.by0 * ((double) NTOTAL);
		fldi.bz[0] = param.bz0 * ((double) NTOTAL);
	}
#endif
}

/** Init the flow arrays... */	
void init_flow(struct Field fldi) {
	int i,n;
	int j,k;
	
	double dummy_var;
	
	DEBUG_START_FUNC;
	// Initialise vectors to 0
	
	if(!param.restart) { //so we don't zero a read-in restart
	  for( n = 0 ; n < fldi.nfield ; n++) {
	    for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
	      fldi.farray[n][i] = 0.0;
	    }
	  }
	}

	if(param.init_large_scale_noise) init_LargeScaleNoise(fldi);

	if(param.init_large_scale_noise_k) init_LargeScaleNoise_k(fldi);
	
	if(param.init_large_scale_2D_noise) init_LargeScale2DNoise(fldi);

	if(param.init_vortex) init_KidaVortex(fldi);

	if(param.init_spatial_structure) init_SpatialStructure(fldi);

	//	if(param.init_Kelvin_Wave) init_Kelvin_Wave(fldi);

	//	if(param.init_Lamb_Vortex) init_Lamb_Vortex(fldi);

	if(param.init_white_noise) init_WhiteNoise(fldi);

	if(param.init_bench) init_Bench(fldi);

	if(param.init_mean_field) init_MeanField(fldi);

	if(param.init_B_field) init_B_field(fldi);
	
	if(param.init_dump) {
		read_dump(fldi, &dummy_var,"init.dmp");
		MPI_Printf("Initial conditions read successfully from the restart dump\n");
	}

#ifdef BOUNDARY_C
	boundary_c(fldi);
#endif
	
       	projector(fldi.vx,fldi.vy,fldi.vz);

#ifdef MHD
	projector(fldi.bx,fldi.by,fldi.bz);
#endif

#ifdef DEBUG
	MPI_Printf("Initflow:\n");
	D_show_all(fldi);
	MPI_Printf("**************************************************************************************\n");
#endif	
	
	DEBUG_END_FUNC;
	
	return;
}
	
	
