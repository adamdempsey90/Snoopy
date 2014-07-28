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


#include <math.h>
#include <complex.h>

#include "common.h"
#include "gfft.h"
#include "debug.h"
#include "forcing.h"

/***************************************************************/
/**
	Compute the right hand side of the INCOMPRESSIBLE dynamical equation
	
	@param dfldo: (output) right hand side of the dynamical equation
	@param fldi: (input) current status of the flow
	@param t: current time of the simulation
	@param tremap: current remap time (only when shear is on)
	@param dt: current timestep size
*/
/***************************************************************/
void timestep( struct Field dfldo,
			   struct Field fldi,
			   const double t,
			   const double tremap,
			   const double dt) {
			   
	int i;
	double complex q0,q1;
	double qr0;
	double S, gamma, epsilon,w;
	double costheta,sintheta;
	int j,k; //AJB

	// This is the timesteping algorithm, solving the physics.

	// Find the shear at time t
#ifdef WITH_SHEAR
#ifdef TIME_DEPENDANT_SHEAR
	S = param.shear * cos(param.omega_shear * t);	// This is the real shear: -dvy/dx
#else
	S = param.shear;
#endif
#endif
#ifdef WITH_ELLIPTICAL_VORTEX //AJB
	gamma = param.gamma;
	epsilon = param.epsilon;
	w = sqrt(gamma*gamma-epsilon*epsilon);
#endif
#ifdef WITH_ROTATION //AJB
	if(param.theta==0.0) {
	  costheta=1.0;
	  sintheta=0.0; //just in case not exact below...
	} else {
	  costheta=cos(param.theta*M_PI/180.0);
	  sintheta=sin(param.theta*M_PI/180.0);
	}
#endif

/* #ifdef ELSASSER_FORMULATION */
/* /\****************************************** */
/* ** ELSASSER variable formulation ********** */
/* ** To be used  */
/* *******************************************\/ */

/* // Solve the MHD equations using Elsasser fields */
/* #ifdef _OPENMP */
/* 	#pragma omp parallel for private(i) schedule(static)	 */
/* #endif */
/* 	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) { */
/* 		w1[i] =  fldi.vx[i]+fldi.bx[i]; */
/* 		w2[i] =  fldi.vy[i]+fldi.by[i]; */
/* 		w3[i] =  fldi.vz[i]+fldi.bz[i]; */
		
/* 		w4[i] =  fldi.vx[i]-fldi.bx[i]; */
/* 		w5[i] =  fldi.vy[i]-fldi.by[i]; */
/* 		w6[i] =  fldi.vz[i]-fldi.bz[i]; */
/* 	} */
	
/* 	// These fields should have no divergence. */
/* 	// When shear is on, however, divergence is conserved up to the timeintegrator precision. */
/* 	// Let's clean it. */
/* 	projector(w1,w2,w3); */
/* 	projector(w4,w5,w6); */

/* 	gfft_c2r_t(w1); */
/* 	gfft_c2r_t(w2); */
/* 	gfft_c2r_t(w3); */
	
/* 	gfft_c2r_t(w4); */
/* 	gfft_c2r_t(w5); */
/* 	gfft_c2r_t(w6); */
	
/* // Compute the Elsasser tensor */

/* #ifdef _OPENMP */
/* 	#pragma omp parallel for private(i) schedule(static)	 */
/* #endif */
/* 	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) { */
/* 		wr7[i]  = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr8[i]  = wr1[i] * wr5[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr9[i]  = wr1[i] * wr6[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr10[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr11[i] = wr2[i] * wr5[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr12[i] = wr2[i] * wr6[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr13[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr14[i] = wr3[i] * wr5[i] / ((double) NTOTAL*NTOTAL); */
/* 		wr15[i] = wr3[i] * wr6[i] / ((double) NTOTAL*NTOTAL); */
/* 	} */
	
/* 	gfft_r2c_t(wr7); */
/* 	gfft_r2c_t(wr8); */
/* 	gfft_r2c_t(wr9); */
/* 	gfft_r2c_t(wr10); */
/* 	gfft_r2c_t(wr11); */
/* 	gfft_r2c_t(wr12); */
/* 	gfft_r2c_t(wr13); */
/* 	gfft_r2c_t(wr14); */
/* 	gfft_r2c_t(wr15); */
	
/* // Compute the volution of the Elssaser fields (u= ik. */
/* #ifdef _OPENMP */
/* 	#pragma omp parallel for private(i) schedule(static)	 */
/* #endif */
/* 	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) { */
/* 		dfldo.vx[i] = - I * 0.5 * mask[i] * ( */
/* 						kxt[i] * ( 2.0 * w7[i]  ) + kyt[i] * ( w8[i] + w10[i]) + kzt[i] * (w9[i]  + w13[i]) ); */
/* 		dfldo.vy[i] = - I * 0.5 * mask[i] * ( */
/* 						kxt[i] * (w10[i] + w8[i]) + kyt[i] * ( 2.0   * w11[i]) + kzt[i] * (w12[i] + w14[i]) ); */
/* 		dfldo.vz[i] = - I * 0.5 * mask[i] * ( */
/* 						kxt[i] * (w13[i] + w9[i]) + kyt[i] * (w14[i] + w12[i]) + kzt[i] * ( 2.0 * w15[i]  ) ); */
		
										
/* 		dfldo.bx[i] = - I * 0.5 * mask[i] * ( */
/* 						                            kyt[i] * ( w8[i] - w10[i]) + kzt[i] * (w9[i]  - w13[i]) ); */
/* 		dfldo.by[i] = - I * 0.5 * mask[i] * ( */
/* 						kxt[i] * (w10[i] - w8[i])                             + kzt[i] * (w12[i] - w14[i]) ); */
/* 		dfldo.bz[i] = - I * 0.5 * mask[i] * ( */
/* 						kxt[i] * (w13[i] - w9[i]) + kyt[i] * (w14[i] - w12[i])  ); */
		

/* 	} */
		
/* // Compute real(U) in case it is used later. */
/* 	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) { */
/* 		wr1[i] = 0.5 * (wr1[i] + wr4[i]); */
/* 		wr2[i] = 0.5 * (wr2[i] + wr5[i]); */
/* 		wr3[i] = 0.5 * (wr3[i] + wr6[i]); */
/* 	} */

/* #else */
/******************************************
** Velocity Self Advection ****************
*******************************************/

		/* Compute the convolution */
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	// These fields should have no divergence.
	// When shear is on, however, divergence is conserved up to the timeintegrator precision.
	// Let's clean it.
	projector(w1,w2,w3);
	
	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
	/* Compute the convolution for the advection process */
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
#ifndef WITH_2D
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
#endif
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
#ifndef WITH_2D
	gfft_r2c_t(wr6);
#endif
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
	  dfldo.vx[i] = - I * mask[i] * ( kxt[i] * w4[i] + kyt[i] * w7[i] + kzt[i] * w8[i] );
	  dfldo.vy[i] = - I * mask[i] * ( kxt[i] * w7[i] + kyt[i] * w5[i] + kzt[i] * w9[i] );
	  dfldo.vz[i] = - I * mask[i] * ( kxt[i] * w8[i] + kyt[i] * w9[i] + kzt[i] * w6[i] );	// since kz=0 in 2D, kz*w6 gives 0, even if w6 is some random array

	} 

	if(param.nonlinearoff) { //AJB turn off nonlinearities
	  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
	    dfldo.vx[i] = 0.0;
	    dfldo.vy[i] = 0.0;
	    dfldo.vz[i] = 0.0;
	  }
	}
	
/* #endif */

/**********************************************
** BOUSSINESQ TERMS (if needed) ***************
***********************************************/

#ifdef BOUSSINESQ
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = fldi.th[i];
	}
	
	gfft_c2r_t(w4);
		
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr5[i] = wr1[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr2[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
#ifndef WITH_2D
		wr7[i] = wr3[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
#endif
#ifdef N2PROFILE
		wr8[i] = N2_profile[i] * wr4[i] / ((double) NTOTAL);
#endif
	}
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
#ifndef WITH_2D
	gfft_r2c_t(wr7);
#endif
#ifdef N2PROFILE
	gfft_r2c_t(wr8);
#endif

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef VERTSTRAT
#ifdef N2PROFILE
		dfldo.vz[i] -= w8[i] * mask[i];
#else
		//dfldo.vz[i] -= param.N2 * fldi.th[i];
                dfldo.vz[i] += fldi.th[i];
#endif
		dfldo.th[i] = - I*mask[i]*(
			kxt[i]*w5[i]+kyt[i]*w6[i]+kzt[i]*w7[i])
			- param.N2*fldi.vz[i];		
		/* dfldo.th[i] = - I * mask[i] * ( */
		/* 	kxt[i] * w5[i] + kyt[i] * w6[i] + kzt[i] * w7[i]) */
		/* 	+ fldi.vz[i]; */
#else //not vertstrat
#ifdef N2PROFILE
		dfldo.vx[i] -= w8[i] * mask[i];
#else
		//	dfldo.vx[i] -= param.N2 * fldi.th[i];
		dfldo.vx[i] += fldi.th[i];
#endif
		dfldo.th[i] = -I*mask[i]*(
		   kxt[i]*w5[i]+kyt[i]*w6[i]+kzt[i]*w7[i])
			- param.N2*fldi.vx[i];		
/* 		dfldo.th[i] = - I*mask[i]*( */
/* 			kxt[i]*w5[i]+kyt[i]*w6[i]+kzt[i]*w7[i]) */
/* 			+ fldi.vx[i]; */
#endif //vertstrat
	}
	
	
#endif

/*********************************************
**** MHD Terms (if needed)   *****************
*********************************************/
#ifdef MHD
#ifndef ELSASSER_FORMULATION		// If Elssaser is on, MHD are already computed...

// Start with the induction equation
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	// These fields should have no divergence.
	// When shear is on, however, divergence is conserved up to the timeintegrator precision.
	// Let's clean it.
	projector(w4,w5,w6);
	
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// (vx,vy,vz) is in w1-w3 and (bx,by,bz) is in (w4-w6). It is now time to compute the emfs in w7-w9...
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = (wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL*NTOTAL);
		wr8[i] = (wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL*NTOTAL);
		wr9[i] = (wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf to add in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.bx[i] = I * mask[i] * (kyt[i] * w9[i] - kzt[i] * w8[i]);
		dfldo.by[i] = I * mask[i] * (kzt[i] * w7[i] - kxt[i]* w9[i]);
		dfldo.bz[i] = I * mask[i] * (kxt[i]* w8[i] - kyt[i] * w7[i]);

#ifdef WITH_ELLIPTICAL_VORTEX //AJB 02/02/12
		dfldo.bx[i] -= (gamma+epsilon)*fldi.by[i];
		dfldo.by[i] -= -(gamma-epsilon)*fldi.bx[i];
#endif

	}


// Let's do the Lorentz Force
// We already have (bx,by,bz) in w4-w6. No need to compute them again...

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
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

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		dfldo.vx[i] += I * mask[i] * (kxt[i] * w1[i] + kyt[i] * w7[i] + kzt[i] * w8[i]);
		dfldo.vy[i] += I * mask[i] * (kxt[i] * w7[i] + kyt[i] * w2[i] + kzt[i] * w9[i]);
		dfldo.vz[i] += I * mask[i] * (kxt[i] * w8[i] + kyt[i] * w9[i] + kzt[i] * w3[i]);
	}
	
#endif
#endif

/************************************
** SOURCE TERMS  ********************
************************************/

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef WITH_ROTATION //AJB modified for tilted rotation in xz-plane
	  dfldo.vx[i]+=2.0*param.omega*fldi.vy[i]*costheta;
	  dfldo.vy[i]-=2.0*param.omega*(fldi.vx[i]*costheta-fldi.vz[i]*sintheta);
	  dfldo.vz[i]-=2.0*param.omega*fldi.vy[i]*sintheta;
#endif
#ifdef WITH_SHEAR
		dfldo.vy[i] += S  * fldi.vx[i];
#ifdef MHD
		dfldo.by[i] -= S * fldi.bx[i];
#endif
#endif
#ifdef WITH_ELLIPTICAL_VORTEX //AJB
		dfldo.vx[i] += (gamma+epsilon)*fldi.vy[i];
		dfldo.vy[i] += -(gamma-epsilon)*fldi.vx[i];
#endif
	}
	
/************************************
** EXPLICIT LINEAR DISSIPATION ******
*************************************/

/* #ifdef WITH_EXPLICIT_DISSIPATION */
/* #ifdef _OPENMP */
/* 	#pragma omp parallel for private(i) schedule(static) */
/* #endif */
/* 	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) { */
/* 		dfldo.vx[i] += - nu * k2t[i] * fldi.vx[i]; */
/* 		dfldo.vy[i] += - nu * k2t[i] * fldi.vy[i]; */
/* 		dfldo.vz[i] += - nu * k2t[i] * fldi.vz[i]; */
		
/* #ifdef MHD */
/* 		dfldo.bx[i] += - eta * k2t[i] * fldi.bx[i]; */
/* 		dfldo.by[i] += - eta * k2t[i] * fldi.by[i]; */
/* 		dfldo.bz[i] += - eta * k2t[i] * fldi.bz[i]; */
/* #endif	// MHD */

/* #ifdef BOUSSINESQ */
/* 		dfldo.th[i] += - nu_th * k2t[i] * fldi.th[i]; */
/* #endif	// BOUSSINESQ */
/* 	} */

/* #endif	// WITH_EXPLICIT_DISSIPATION */
	
/************************************
** PRESSURE TERMS *******************
************************************/

#ifdef _OPENMP
	#pragma omp parallel for private(i,q0,q1) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			
#ifdef WITH_SHEAR
	  q0= S * ky[i] * fldi.vx[i] + kxt[i] * dfldo.vx[i] + ky[i] * dfldo.vy[i] + kzt[i] * dfldo.vz[i];
#ifdef WITH_ELLIPTICAL_VORTEX //AJB
	    q0 = (gamma+epsilon)*kxt[i]*fldi.vy[i]-(gamma-epsilon)*kyt[i]*fldi.vx[i] + kxt[i]*dfldo.vx[i] + kyt[i]*dfldo.vy[i] + kzt[i]*dfldo.vz[i];	  
#endif
#else
	  q0= kxt[i] * dfldo.vx[i] + kyt[i] * dfldo.vy[i] + kzt[i] * dfldo.vz[i];
#endif
/* po would contain the pressure field
		if(po != NULL) {
			po[i] = - I * ik2t[i] * q0;	// Save the pressure field (if needed)
		}
*/
		dfldo.vx[i] += -kxt[i]* q0 * ik2t[i];
		dfldo.vy[i] += -kyt[i]* q0 * ik2t[i];
		dfldo.vz[i] += -kzt[i]* q0 * ik2t[i];
	}

	return;
}

/***************************************************************/
/**
	Implicit steps of the integrator (essentially linear diffusion terms)
	
	This is an implicit model: fldi is modified by this routine
	
	@param fldi: (input and output) current status of the flow
	@param t: current time of the simulation
	@param dt: current timestep size
*/
/***************************************************************/
void implicitstep(
			   struct Field fldi,
			   const double t,
			   const double dt ) {
			   
  double q0;
  double flux,tstartrelax,fluxtol,tau,avgfluxnow,vol;
  int i;

#ifndef WITH_EXPLICIT_DISSIPATION
#ifdef _OPENMP
	#pragma omp parallel for private(i,q0) schedule(static)
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {

	  if(param.hyperdiffusion) {//nabla^4
	    q0 = exp( - nu * dt* pow(k2t[i],2.0) );
	  } else {//nabla^2
	    q0 = exp( - nu * dt* k2t[i] );
	  }

	  fldi.vx[i] = fldi.vx[i] * q0;
	  fldi.vy[i] = fldi.vy[i] * q0;
	  fldi.vz[i] = fldi.vz[i] * q0;
		
#ifdef BOUSSINESQ
	  if(param.hyperdiffusion) {//nabla^4
	    q0 = exp( - nu_th * dt* pow(k2t[i],2.0) ); //AJB 19/04/13
	    //q0 = exp( - nu_th * dt* k2t[i] );
	  } else {//nabla^2
	    q0 = exp( - nu_th * dt* k2t[i] );
	  }
	  fldi.th[i] = fldi.th[i] * q0;
	  if(param.heating){//AJB heating/cooling using integratig factor
	    if(!param.hyperdiffusion) {
	      fldi.th[i] += w15[i]*(1.0-q0)*ik2t[i]/nu_th; 
	    } else {//hyperdiffusion
	      fldi.th[i] += w15[i]*(1.0-q0)*pow(ik2t[i],2.0)/nu_th;
	    }
	  } //AJB
#endif
#ifdef MHD
	  if(param.hyperdiffusion) {//nabla^4
	    q0 = exp( - eta * dt* pow(k2t[i],2.0) );
	  } else {//nabla^2
	    q0 = exp( - eta * dt* k2t[i] );
	  }
	  fldi.bx[i] = fldi.bx[i] * q0;
	  fldi.by[i] = fldi.by[i] * q0;
	  fldi.bz[i] = fldi.bz[i] * q0;
#endif
	}
#endif	// WITH_EXPLICIT_DISSIPATION

#ifdef FORCING
	forcing(fldi, dt);
#endif
	if(param.fluxrelax){ 	//AJB 04/07/13
	  tstartrelax=2.0; fluxtol=0.02; tau=500.0; //may need tweaking...
	  if(t>tstartrelax){//allow convection to settle a little first
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(static)
#endif
	    for(i=0;i<NTOTAL_COMPLEX;i++){
	      w1[i]=fldi.th[i]; //b
	      w2[i]=fldi.vz[i]; //uz
	    }
	    gfft_c2r_t(w1); gfft_c2r_t(w2);
	    flux=0.0;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(static) reduction(+:flux)
#endif
	    for(i=0;i<2*NTOTAL_COMPLEX;i++){
	      wr1[i]/=((double) NTOTAL); wr2[i]/=((double) NTOTAL);
	      flux+=(wr1[i]*wr2[i])/((double) NTOTAL);
	    }
	    reduce(&flux,1); //mean convective flux
	    vol=(param.lz/2.0); //lx,ly accounted for in other terms...
	    param.flux_runningavg+=(flux-param.N2*nu_th)*dt*vol; //total flux
	    avgfluxnow=param.flux_runningavg/(t-tstartrelax);
	    //set flux=1 over timescale tau, switching off when within few %
	    if(fabs(avgfluxnow+(-(param.lz/2.0)+(1.0+nu_th*param.alpha)*2.0*param.zc)/vol)>fluxtol){//relax until reached tolerance
	      if(!param.hyperdiffusion){
		param.N2+=(((flux+(-(param.lz/2.0)+(1.0+nu_th*param.alpha)*2.0*param.zc)/vol)/nu_th)-param.N2)*dt/tau;
	      } else {//hyperdiffusion means 1/nu_th is large...
		param.N2+=((flux+(-(param.lz/2.0)+2.0*param.zc)/vol)/0.001)*dt/tau;
	      } 
	    }
//	    MPI_Printf("t,N2,F,uzb,N2target,F_runningavg,F_target: %f,%f,%f,%f,%f,%f,%f \n",t,param.N2,flux-param.N2*nu_th,flux,(flux+(-(param.lz/2.0)+(1.0+nu_th*param.alpha)*2.0*param.zc)/vol)/nu_th,avgfluxnow/vol,((param.lz/2.0)-(1.0+nu_th*param.alpha)*2.0*param.zc)/vol);
	  } 
	} //end AJB
	return;
}
