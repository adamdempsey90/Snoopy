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
#include "rk45.h"
#ifdef MPI_SUPPORT
#include "transpose.h"
#endif
#include "debug.h"
#ifdef WITH_SHEAR

double time_shift(double t) {
	double tremap;
#ifdef TIME_DEPENDANT_SHEAR
	tremap = sin(param.omega_shear * t) / param.omega_shear;	// This is formally the angular displacement divded by param.shear= int dt S(t) /<S>
#else
	tremap = fmod(t + param.ly / (2.0 * param.shear * param.lx) , param.ly / (param.shear * param.lx)) - param.ly / (2.0 * param.shear * param.lx);
#endif
	return(tremap);
}

void remap(double complex qi[]) {
	int i, j, k;
	int nx, ny, nxtarget;
	
	DEBUG_START_FUNC;
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i]=0.0;
	}
	
#ifdef MPI_SUPPORT
// We have to transpose the array to get the remap properly
	transpose_complex_XY(qi,qi);
	
	for( i = 0; i < NX_COMPLEX; i++) {
		nx = fmod( i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 ;
		for( j = 0; j < NY_COMPLEX/NPROC; j++) {
			ny = fmod( j + rank * NY_COMPLEX / NPROC + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 ;
			
			nxtarget = nx + ny;		// We have a negative shear, hence nx plus ny
			
			if( (nxtarget > -NX_COMPLEX / 2) & (nxtarget < NX_COMPLEX/2)) {
			
				if ( nxtarget <0 ) nxtarget = nxtarget + NX_COMPLEX;
			
				for( k = 0; k < NZ_COMPLEX; k++) {
					w1[k + NZ_COMPLEX * nxtarget + NZ_COMPLEX * NX_COMPLEX * j] = qi[ k + i * NZ_COMPLEX + j * NZ_COMPLEX * NX_COMPLEX];
				}
			}
		}
	}
	
	// transpose back
	transpose_complex_YX(w1,w1);

#else
	for( i = 0; i < NX_COMPLEX; i++) {
		nx = fmod( i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 ;
		for( j = 0; j < NY_COMPLEX; j++) {
			ny = fmod( j + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 ;
			
			nxtarget = nx + ny;		// We have a negative shear, hence nx plus ny
			
			if( (nxtarget > -NX_COMPLEX / 2) & (nxtarget < NX_COMPLEX/2)) {
			
				if ( nxtarget <0 ) nxtarget = nxtarget + NX_COMPLEX;
			
				for( k = 0; k < NZ_COMPLEX; k++) {
					w1[k + NZ_COMPLEX * j + NZ_COMPLEX * NY_COMPLEX * nxtarget] = qi[ IDX3D ];
				
				}
			}
		}
	}
#endif

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		qi[i] = w1[i] * mask[i];
	}

	DEBUG_END_FUNC;
	
	return;
}

void kvolve(const double tremap) {
	int i, j, k;
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k) schedule(static)	
#endif	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kxt[ IDX3D ] = kx[ IDX3D ] + tremap * param.shear * ky[ IDX3D ];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
						       ky[IDX3D] * ky[IDX3D]+
							   kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}

	return;
}
#endif
//AJB
#ifdef WITH_ELLIPTICAL_VORTEX
void kvolvellipse(const double t) {
  int indx;
	double E, invE, w;
	E = sqrt((param.gamma+param.epsilon)/(param.gamma-param.epsilon));
	invE = 1.0/E;
	w = sqrt(param.gamma*param.gamma-param.epsilon*param.epsilon);
	if(param.gamma == 0.0) {//just to make sure...
	  E = 0.0; w = 0.0; invE = 0.0;
	}	
#ifdef _OPENMP
#pragma omp parallel for private(indx) schedule(static)
#endif
	for (indx=0; indx < NTOTAL_COMPLEX; indx++) {
	  kxt[indx] = kx[indx]*cos(w*t) - invE*ky[indx]*sin(w*t);
	  kyt[indx] = E*kx[indx]*sin(w*t) + ky[indx]*cos(w*t);
	  k2t[ indx ] = kxt[indx] * kxt[indx] +
	    kyt[indx] * kyt[indx]+
	    kz[indx] * kz[indx];
	  if ( k2t[indx] == 0.0 ) ik2t[indx] = 1.0;
	  else	ik2t[indx] = 1.0 / k2t[indx];
	}	
  return;
}
void knumvolvellipse(const double t, const double dt, const int kflag) {
  int i, j, k, indx, flag;
  double y[9], yp[9];
  double relerr, abserr;

  relerr = 1e-20; abserr = 1e-20;
  for(i=0;i<9;i++){ y[i] = kbasis[i]; yp[i] = 0.0;}
  //MPI_Printf("%15.15e \t %15.15e \t %15.15e \n",y[0],y[1],y[2]);
  flag = r8_rkf45 ( f , 9, y, yp, &t, t+dt, &relerr, abserr, 1);
  //MPI_Printf("%15.15e \t %15.15e \t %15.15e \n",y[0],y[1],y[2]);
  if(kflag == 1) { //if END of full RK timestep store new value (NOT otherwise!!!!)
    for(i=0;i<9;i++){ kbasis[i] = y[i];  } 
  }

#ifdef _OPENMP
#pragma omp parallel for private(indx) schedule(static)	
#endif	
  for (indx=0; indx < NTOTAL_COMPLEX; indx++) {
    kxt[indx] = kx[indx]*y[0] + ky[indx]*y[3] + kz[indx]*y[6];
    kyt[indx] = kx[indx]*y[1] + ky[indx]*y[4] + kz[indx]*y[7];
    kzt[indx] = kx[indx]*y[2] + ky[indx]*y[5] + kz[indx]*y[8];
    
    k2t[ indx ] = kxt[indx] * kxt[indx] +
      kyt[indx] * kyt[indx] +
      kzt[indx] * kzt[indx];
    
    if ( k2t[indx] == 0.0 ) ik2t[indx] = 1.0;
    else	ik2t[indx] = 1.0 / k2t[indx];
  }
  
  return;
}
#endif
/******************************************************************************/
//Evaluate derivates i.e. RHS of dot{k} = -A^{T} k
void f ( double t, double y[], double yp[] ) {
  yp[0] = (param.epsilon-param.gamma)*y[1]; //k1[1]
  yp[1] = (param.gamma+param.epsilon)*y[0]; //k1[2]
  yp[2] = 0.0; //k1[3]
  yp[3] = (param.epsilon-param.gamma)*y[4]; //k2[1]
  yp[4] = (param.gamma+param.epsilon)*y[3]; //k2[2]
  yp[5] = 0.0; //k2[3]
  yp[6] = 0.0; //k3[1]
  yp[7] = 0.0; //k3[2]
  yp[8] = 0.0; //k3[3]
  return;
}
