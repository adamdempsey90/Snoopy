#include <stdlib.h>
#include <string.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../snoopy.h" //AJB

/***********************************************************/
/** 
	find the maximum of a real array of size 2*NTOTAL_COMPLEX
	@param wri array in which we want to know the maximum
*/
/***********************************************************/

double find_max(double *wri) {
	double q0;
	int i,j,k,idx;
	q0=wri[0];
	
	for(i=0 ; i < NX/NPROC ; i++) {
		for(j=0 ; j < NY ; j++) {
			for(k=0 ; k < NZ ; k++) {
#ifdef WITH_2D
				idx = j + i * (NY+2);
#else
				idx = k + j * (NZ + 2) + i * NY * (NZ + 2);
#endif
				if(q0<wri[idx]) q0=wri[idx];
			}
		}
	}
	
	return(q0);
}

/***********************************************************/
/** 
	find the minimum of a real array of size 2*NTOTAL_COMPLEX
	@param wri array in which we want to know the minimum
*/
/***********************************************************/

double find_min(double *wri) {
	double q0;
	int i,j,k,idx;
	q0=wri[0];
	
	for(i=0 ; i < NX/NPROC ; i++) {
		for(j=0 ; j < NY ; j++) {
			for(k=0 ; k < NZ ; k++) {
#ifdef WITH_2D
				idx = j + i * (NY+2);
#else
				idx = k + j * (NZ + 2) + i * NY * (NZ + 2);
#endif
				if(q0>wri[idx]) q0=wri[idx];
			}
		}
	}
	
	return(q0);
}

/***********************************************************/
/** 
	compute the correlation between 2 fields
*/
/***********************************************************/
double compute_2correlation(double *wri1, double *wri2) {
	double q0;
	int i,j,k;
	q0=0.0;
	
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		q0 += wri1[i] * wri2[i] / ((double) NTOTAL);
	}
	
	return(q0);
}

/***********************************************************/
/** 
	compute the correlation between 3 fields
*/
/***********************************************************/
double compute_3correlation(double *wri1, double *wri2, double *wri3) {
	double q0;
	int i,j,k;
	q0=0.0;
	
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		q0 += wri1[i] * wri2[i] * wri3[i] / ((double) NTOTAL);
	}
	
	return(q0);
}

/***********************************************************/
/** 
	Write statistical quantities using text format in the file
	timevar. 
	@param fldi Field structure from which the statistical quantities are derived.
	@param t Current time of the simulation
*/
/***********************************************************/

void output_timevar(const struct Field fldi,
					const double t) {
					
	FILE *ht;
	double output_var;
	static int warning_flag = 0;

	int i,j,k,tot,var;
        double z,temp1,temp2;
	
	DEBUG_START_FUNC;
	
	// Open the timevar file
	if(rank==0) {
	  ht=fopen(OUTPUT_TIMEVAR_FILENAME,"a");
	}
	
	for(i=0;i<NTOTAL_COMPLEX;i++){
	  w1[i]=fldi.vx[i]; //ux
	  w2[i]=fldi.vy[i]; //uy
	  w3[i]=fldi.vz[i]; //uz
	}
	gfft_c2r_t(w1); gfft_c2r_t(w2); gfft_c2r_t(w3);
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
	  wr1[i] = wr1[i] / ((double) NTOTAL ); //ux
	  wr2[i] = wr2[i] / ((double) NTOTAL ); //uy
	  wr3[i] = wr3[i] / ((double) NTOTAL ); //uz
	}

	// loop on all the requested variables
	for( var = 0 ; var < param.timevar_vars.length ; var++ ) {
	        if(!strcmp(param.timevar_vars.name[var],"t")) {
			// current time
		  output_var=t;
		}
		else if(!strcmp(param.timevar_vars.name[var],"ev")) {
			// kinetic energy
		  output_var = energy(fldi.vx) + energy(fldi.vy) + energy(fldi.vz);
		  reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"ed")) {
			// kinetic energy dissipation due to viscosity. Factor of 2 because KE is 1/2 u^2
		  output_var = (2.0*nu)*(dissipation_energy(fldi.vx)+dissipation_energy(fldi.vy)+dissipation_energy(fldi.vz));
		  reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"vxmax")) {
			// maximum of vx component
			output_var=find_max(wr1);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[var],"vymax")) {
			// maximum of vy component
			output_var=find_max(wr2);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[var],"vzmax")) {
			// maximum of vz component
			output_var=find_max(wr3);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[var],"vxmin")) {
			// minimum of vx component
			output_var=find_min(wr1);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[var],"vymin")) {
			// minimum of vy component
			output_var=find_min(wr2);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[var],"vzmin")) {
			// minimum of vz component
			output_var=find_min(wr3);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[var],"vxvy")) {
			// incompressible reynolds stress
			output_var=compute_2correlation(wr1,wr2);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"uAu")) { //AJB
			// total incompressible reynolds stress (ONLY FOR ELLIPTICAL VORTEX PROBLEM IN ROTATING FRAME)
		        // same as vxvy in non-rotating frame
		  output_var=sin(2.0*param.gamma*t)*2.0*(energy(fldi.vx)-energy(fldi.vy)) + 2.0*cos(2.0*param.gamma*t)*compute_2correlation(wr1,wr2);
		  reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"hv")) {
			// kinetic helicity
			// Compute vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * ik2t[j] * (kyt[j] * fldi.vz[j] - kzt[j] * fldi.vy[j] );
				w11[j] = I * ik2t[j] * (kzt[j] * fldi.vx[j] - kxt[j] * fldi.vz[j] );
				w12[j] = I * ik2t[j] * (kxt[j] * fldi.vy[j] - kyt[j] * fldi.vx[j] );
			}
	
			gfft_c2r(w10);
			gfft_c2r(w11);
			gfft_c2r(w12);
	
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr10[j] = wr10[j] / ((double) NTOTAL );
				wr11[j] = wr11[j] / ((double) NTOTAL );
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}

			output_var=compute_2correlation(wr10,wr1)+compute_2correlation(wr11,wr2)+compute_2correlation(wr12,wr3);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"w2")) {
			// total enstrophy
			// Compute vorticity
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * (kyt[j] * fldi.vz[j] - kzt[j] * fldi.vy[j]);
				w11[j] = I * (kzt[j] * fldi.vx[j] - kxt[j] * fldi.vz[j]);
				w12[j] = I * (kxt[j] * fldi.vy[j] - kyt[j] * fldi.vx[j]);
			}
			output_var=energy(w10)+energy(w11)+energy(w12);
			reduce(&output_var,1);
		}
#ifdef MHD
		else if(!strcmp(param.timevar_vars.name[var],"edm")) {
			// kinetic energy dissipation due to magnetic diffusion. Factor of 2 because KE is 1/2 u^2
		  output_var = (2.0*eta)*(dissipation_energy(fldi.bx)+dissipation_energy(fldi.by)+dissipation_energy(fldi.bz));
		  reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"em")) {
			// magnetic energy
			output_var = energy(fldi.bx) + energy(fldi.by)+energy(fldi.bz);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"bxmax")) {
			// maximum of bx component
			output_var=find_max(wr5);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[var],"bymax")) {
			// maximum of by component
			output_var=find_max(wr6);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[var],"bzmax")) {
			// maximum of bz component
			output_var=find_max(wr7);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[var],"bxmin")) {
			// minimum of bx component
			output_var=find_min(wr5);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[var],"bymin")) {
			// minimum of by component
			output_var=find_min(wr6);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[var],"bzmin")) {
			// minimum of bz component
			output_var=find_min(wr7);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[var],"bxby")) {
			// incompressible maxwell stress
			output_var=compute_2correlation(wr5,wr6);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"hc")) {
			// cross helicity (u.b)
			output_var=compute_2correlation(wr1,wr5)+compute_2correlation(wr2,wr6)+compute_2correlation(wr3,wr7);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"hm")) {
			// magnetic helicity
			// Compute vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * ik2t[j] * (kyt[j] * fldi.bz[j] - kzt[j] * fldi.by[j] );
				w11[j] = I * ik2t[j] * (kzt[j] * fldi.bx[j] - kxt[j]* fldi.bz[j] );
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - kyt[j] * fldi.bx[j] );
			}
	
			gfft_c2r(w10);
			gfft_c2r(w11);
			gfft_c2r(w12);
	
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr10[j] = wr10[j] / ((double) NTOTAL );
				wr11[j] = wr11[j] / ((double) NTOTAL );
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}

			output_var=compute_2correlation(wr10,wr5)+compute_2correlation(wr11,wr6)+compute_2correlation(wr12,wr7);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"j2")) {
			// total current rms
			// Compute current
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * (kyt[j]  * fldi.bz[j] - kzt[j]  * fldi.by[j]);
				w11[j] = I * (kzt[j]  * fldi.bx[j] - kxt[j] * fldi.bz[j]);
				w12[j] = I * (kxt[j] * fldi.by[j] - kyt[j]  * fldi.bx[j]);
			}
			output_var=energy(w10)+energy(w11)+energy(w12);
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[var],"az2")) {
			// Square of the vertical component of the vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - kyt[j] * fldi.bx[j] );
			}
			
			output_var=energy(w12);
			
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[var],"vxaz")) {
			// Source term in the streamfunction equation
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - kyt[j] * fldi.bx[j] );
			}
			
			gfft_c2r(w12);
			
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}
			
			output_var=compute_2correlation(wr12,wr1);
			reduce(&output_var,1);
			
		}
#endif
#ifdef BOUSSINESQ
		else if(!strcmp(param.timevar_vars.name[var],"et")) { //AJB make sure this is called first
		  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		    w4[i] = fldi.th[i];
		  }
		  gfft_c2r_t(w4);
		  for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		    wr4[i] = wr4[i] / ((double) NTOTAL );
		  }
			// thermal energy
			output_var = energy(fldi.th);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"thmax")) {
			// maximum of th
			output_var=find_max(wr4);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[var],"thmin")) {
			// minimum of th
			output_var=find_min(wr4);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[var],"thvx")) {
			// turbulent heat flux in x
			output_var=compute_2correlation(wr4,wr1);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"thvz")) {
		  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		    w4[i] = fldi.th[i];
		  }
		  gfft_c2r_t(w4);
		  for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
                    wr4[i] = wr4[i] / ((double) NTOTAL );
		  }
		  // turbulent heat flux in z
		  output_var=compute_2correlation(wr4,wr3);
		  reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"delT")) {
		  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		    w4[i]=fldi.th[i];}
		  gfft_c2r_t(w4);
		  for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		    wr4[i]/=((double) NTOTAL );}
		  temp1=0.0; temp2=0.0;
		  for(i = 0 ; i < NX/NPROC ; i++) {
		    for(j = 0 ; j < NY ; j++) {
		      k=0; //z=0
		      temp1+=wr4[k+(NZ+2)*j+(NZ+2)*NY*i]/=((double) NX*NY);
		      k=NZ/2; //z=L
		      temp2+=wr4[k+(NZ+2)*j+(NZ+2)*NY*i]/=((double) NX*NY);
		    } }
		  //mean temperature drop across box
		  reduce(&temp1,1);
		  reduce(&temp2,1);
		  output_var=temp2-temp1;
		  reduce(&output_var,1); 
		  output_var/=((double) NPROC);
		}
	        else if(!strcmp(param.timevar_vars.name[var],"uzrms")) {
			// turbulent heat flux in z
		        output_var = energy(fldi.vz);
			reduce(&output_var,1);
			output_var = sqrt(2.0*output_var);
			reduce(&output_var,1); 
			output_var/=((double) NPROC);
		}
		else if(!strcmp(param.timevar_vars.name[var],"uzrms")) {
			// turbulent heat flux in z
		        output_var = energy(fldi.vz);
			reduce(&output_var,1);
			output_var = sqrt(2.0*output_var);
			reduce(&output_var,1); 
			output_var/=((double) NPROC);
		}
		else if(!strcmp(param.timevar_vars.name[var],"thvzmid")) {
		  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		    w5[i] = fldi.th[i];
		  }
		  gfft_c2r_t(w5);
		  for(i = 0 ; i < NX/NPROC ; i++) {
		    for(j = 0 ; j < NY ; j++) {
		      k=NZ/2.0; 
		      wr5[k+(NZ+2)*j+(NZ+2)*NY*i]/=((double) NTOTAL);
		      output_var+=wr3[k+(NZ+2)*j+(NZ+2)*NY*i]*wr5[k+(NZ+2)*j+(NZ+2)*NY*i]/((double) NX*NY);
		    }
		  }
		  reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[var],"uzmid")) {
		  for(i = 0 ; i < NX/NPROC ; i++) {
		    for(j = 0 ; j < NY ; j++) {
		      k=NZ/4.0;
		      output_var+=wr3[k+(NZ+2)*j+(NZ+2)*NY*i]*wr3[k+(NZ+2)*j+(NZ+2)*NY*i]/((double) NX*NY);
		    }
		  }
		  reduce(&output_var,1);
		  output_var=sqrt(output_var);
		  reduce(&output_var,1); 
		  output_var/=((double) NPROC);
		}
		else if(!strcmp(param.timevar_vars.name[var],"N2mid")) {
                  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		    w5[i] = I*kzt[i]*fldi.th[i];
		  }
		  gfft_c2r_t(w5);
		  for(i = 0 ; i < NX/NPROC ; i++) {
		    for(j = 0 ; j < NY ; j++) {
		      k=NZ/4.0; //z=Lz/2 in physical space
		      wr5[k+(NZ+2)*j+(NZ+2)*NY*i]/=((double) NTOTAL);
		      output_var+=wr5[k+(NZ+2)*j+(NZ+2)*NY*i]/((double) NX*NY);
		    }
		  }
		  reduce(&output_var,1);
		}
	        else if(!strcmp(param.timevar_vars.name[var],"thtotvz")) {
                  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		    w5[i] = fldi.th[i];
		  }
		  gfft_c2r_t(w5);
		  for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
                    wr5[i] = wr5[i] / ((double) NTOTAL );
		  }
		  // turbulent heat flux in z
		  output_var=compute_2correlation(wr5,wr3);
		  reduce(&output_var,1);
		  output_var-=param.N2*nu_th;
		  reduce(&output_var,1); 
		  output_var/=((double) NPROC);
		}
		else if(!strcmp(param.timevar_vars.name[var],"khpeakE")) {
                  for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		    w14[i]=pow(kxt[i]*kxt[i]+kyt[i]*kyt[i],0.25)*fldi.vx[i];
		    w16[i]=pow(kxt[i]*kxt[i]+kyt[i]*kyt[i],0.25)*fldi.vy[i];
		    w17[i]=pow(kxt[i]*kxt[i]+kyt[i]*kyt[i],0.25)*fldi.vz[i];
		  }
		  temp1=energy(w14)+energy(w16)+energy(w17);
		  reduce(&temp1,1); 
		  temp2=energy(fldi.vx)+energy(fldi.vy)+energy(fldi.vz);
		  reduce(&temp2,1);
		  output_var=temp1/temp2;
		  reduce(&output_var,1); output_var/=((double) NPROC);
		}
#endif
		else {
			if(!warning_flag) {
				ERROR_HANDLER(ERROR_WARNING,"Unable to produce the requested data\n");
				MPI_Printf("Timevar output string ''%s'' is unknown\n",param.timevar_vars.name[var]);
				warning_flag=1;
			}
			output_var=0.0;
		}
		
		// Write the variable in the timevar file
		if(rank==0) {
			fprintf(ht,"%08e\t",output_var);
			if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing timevar file");
		}
		
#ifdef MPI_SUPPORT
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
		
	if(rank==0) {
		fprintf(ht,"\n");
		fclose(ht);
	}

	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	Remove the timevar file (if exists) to start from a fresh one.
*/
/**************************************************************************************/
void init_timevar() {
	FILE* ht;
	int i;
	DEBUG_START_FUNC;
	
	if(rank==0) {

	  //ht=fopen("data/e0p2i1j1k1/timevar","w");
	  ht=fopen(OUTPUT_TIMEVAR_FILENAME,"w");

		// print a line with the fields we are going to write
		
		for( i = 0 ; i < param.timevar_vars.length ; i++) {
			fprintf(ht,"%s\t\t",param.timevar_vars.name[i]);
		}
		fprintf(ht,"\n");
		
		fclose(ht);
	}

	DEBUG_END_FUNC;
	return;
}


