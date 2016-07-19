typedef struct {
  
  double FxdPars[30];
  
  
} ParamStruct;

int Debug=0;

#include <time.h>
#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>	// GSL_SUCCEnn ...
#include <gsl/gsl_odeiv.h>	// ODE solver
#include <gsl/gsl_odeiv2.h>	// ODE solver

#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "math.h"
#include "stdio.h"

//Custom Header, for integration routine:
#include "JOE_fast_odesGD.h"

//Custom functions:
#include "Nrutil.c"
#include "Integration_Func.c"
#include "Dispersal_Func.c"
//#include "JOE_DistdRK45.c"

long seed;
gsl_rng *r2;  //This has to be global, to ensure that the generator doesn't start over again.


int main(void){
  
  //Basic parameters, for loops, file printing, random numbers, etc.
  FILE *fp, *fp1;
  int i;
  
  // seed
  long seed;
  srand((unsigned) time(NULL));
  seed = time(NULL)*(int)getpid();
  
  gsl_rng_env_setup ();
  gsl_rng_default_seed = seed;
  const gsl_rng_type *T2;
  T2 = gsl_rng_default;
  r2 = gsl_rng_alloc(T2);
  
  /**********************************/
  /* Within-season model parameters */
  /**********************************/
  
  // FxdPars[0]: NULL
  // FxdPars[1]: 1/C^2 (k)
  // FxdPars[2]: nu bar trans
  // FxdPars[3]: mu (decay)
  // FxdPars[4]: sigma - actually var
  // FxdPars[5]: ratio
  // FxdPars[6]: delta
  // FxdPars[7]: m, #E stages
  // FxdPars[8]: feed (feeding rate)
  // FxdPars[9]: bgrow (Biomass growth rate)
  // FxdPars[10]: S0 (initial suscept host density)
  // FxdPars[11]: Z0 (initial virus/cadaver density)
  // FxdPars[12]: Biom0 (initial Biomass denisty)
  // FxdPars[13]: Fol0 (initial foliage density)
  // FxdPars[14]: 
  // FxdPars[15]: Storage for ODE (S0)
  
  
  double FxdPars[30];
  // Scan in the SEIR model params from a file
  int PCount = 7;
  fp = fopen("Params.dat", "r");
  for(i=1;i<=PCount;i++){
    fscanf(fp, "%lf\t", &FxdPars[i]);
  }
  fclose(fp);
  
  int MaxDays = 60; // Number of days in epizootic 
  int reps = 1; // Number of stochastic realizations per plot
  
  
  /*Biomass and Foliage Within-season Growth*/

  double feed = 0.03; //feeding rate of larvae on foliage
  double bgrow = 0.00005; // woody biomass growth rate
  
  FxdPars[8] = feed;
  FxdPars[9] = bgrow;
  
  
  /* INITIAL CONDITIONS */
  /* These will be log-normal, for long tails */
  
  double lmean_N0, lmean_Z0, lmean_B0; // means on log scale
  double lsd_N0, lsd_Z0, lsd_B0; // sd of the log-normal dist
  
  // For now, the specific values are more or less arbitrary
  lmean_N0 = log(5);
  lmean_Z0 = log(0.05);
  lmean_B0 = log(1500);
  
  lsd_N0 = 0.75;
  lsd_Z0 = 0.75;
  lsd_B0 = 0.25;
  
  
  /********************************/  /********************************/
  /********************************/  /********************************/
  

  /************************************************/
  /*Between-year Parameters (Host, Virus and Tree)*/
  /************************************************/
  
  
  // Longitudinal simulation:
  int MaxYear = 4; // How many years
  
  //Host dynamics:
  double reprod = exp(0.962); //0.962=r from Mason1974 lambda = exp(r)
  double winter = 8.0;
  double gamma = 0.1; // long-term persistence of Z, range (0,1)
  double a_pred = 0; //0 means no predator model
  double b_pred = 0.16;
  
  double disperse = 0.1; //average percent of individuals that disperse
  
  //Foliage model params:
  double fgrow_winter = 0.05; //As a fraction of total branch biomass
  
  //Biomass model params:
  double bgrow_winter = 2;
  double a_biom = 1.15; //a and c determine loss of biomass due to defoliation
  double c_biom = 0.2;
  
  
  /********************************/  /********************************/
  /********************************/  /********************************/
  
  
  /**********************************************************/
  /* Set up the Square Lattice                              */
  /* Edge cells have a different number of neighbors        */
  /**********************************************************/
  
  double Lat_size = 8; //Number of rows and columns 
  double Cells = pow(Lat_size, 2); //Total number of cells in lattice

  //Storage:
  float *N, *Z, *B, *F;
  N = vector(1, Cells); 
  Z = vector(1, Cells); 
  B = vector(1, Cells); 
  F = vector(1, Cells);

  
  //Assign initial values for each cell:
  
  int c;
  
  for(c=1;c<=Cells;c++){
    
    N[c] = exp( gsl_ran_gaussian(r2, lsd_N0) + lmean_N0 );
    Z[c] = exp( gsl_ran_gaussian(r2, lsd_Z0) + lmean_Z0 );
    B[c] = exp( gsl_ran_gaussian(r2, lsd_B0) + lmean_B0 );
    // Foliage is going to be a fraction of the available biomass:
    F[c] = gsl_ran_flat(r2, 0.3, 0.8) * B[c];
    
  }
  
  
  
  /********************************/  /********************************/
  /********************************/  /********************************/
  
  
  /**************************************************************/
  /* Implement the Longitudinal model                           */
  /* Relies on outside functions for integration and dispersal  */
  /**************************************************************/
  
  
  int y;
  int year, cell;
  float *Int_Out; //Output of integration function
  
  for(y=1;y<=MaxYear;y++){
    
    float *oldN, *oldZ, *oldB, *oldF;
    float *Z_begin, *Percent_Defol;
    oldN = vector(1, Cells);
    oldZ = vector(1, Cells); 
    oldB = vector(1, Cells);
    oldF = vector(1, Cells);
    Z_begin = vector(1, Cells); 
    Percent_Defol = vector(1, Cells);
    
    if(Debug == 1){
      printf("Start SEIR integration for all cells\n");
    }
    //Run the within-season model for each cell:
    for(c=1;c<=Cells;c++){
      
      year = y;
      cell = c;
      
      Z_begin[c] = Z[c]; // this is for the overwinter model for Z
      
      //Run the within-season model:
      Int_Out = Integration_Func(FxdPars, N, Z, B, F, year, cell, MaxDays, r2);
      
      // Record the new, old values:
      oldN[c] = *(Int_Out+0);
      oldZ[c] = *(Int_Out+1);
      oldB[c] = *(Int_Out+2);
      oldF[c] = *(Int_Out+3);
      Percent_Defol[c] = *(Int_Out+4);
      
    } //End cell SEIR
    
    if(Debug == 1){
      printf("End SEIR integration for all cells\n");
    }
    
    /******************
     * Host dispersal
     * Tree re-seeding 
     ******************/
    if(Debug == 1){
      printf("Starting dispersal function \n");
    }
    
    Dispersal_Func(oldN,oldZ,oldB,
                   Lat_size,Cells,disperse,r2, y);
    
    if(Debug == 1){
      printf("End dispersal function \n");
    }
    
    /*************************
     * OVERWINTER REPRODUCTION
     *************************/
    if(Debug == 1){
      printf("Start overwintering \n");
    }
    
    for(c=1;c<=Cells;c++){
      
      //Print the beginning and ending values:
      fp1 = fopen("Output.dat","a");
      fprintf(fp1, "%d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
              y, c,
              N[c], Z[c], B[c], F[c], 
              oldN[c], oldZ[c], oldB[c], oldF[c],
              Percent_Defol[c]
              );
      fclose(fp1);
      
      /*************************
      * OVERWINTER REPRODUCTION
      *************************/
      
      //Calculate starting conditions for N w/ predator model:
      N[c] = oldN[c] * reprod; //* (1 - ((a_pred*b_pred*N_start)/(pow(b_pred,2)+pow(N_start,2))));
      
      //Calculate starting conditions for Z:
      Z[c] = oldZ[c] * winter + gamma * Z_begin[c]; //This Z_begin is from the previous generation, so Z_n
      
      //Calculate starting conditions for Biomass:
      B[c] = oldB[c] * bgrow_winter * (1 - (a_biom * Percent_Defol[c] / (c_biom + Percent_Defol[c])));
      
      //Calculate starting conditions for Foliage:
      F[c] = oldF[c] + oldB[c] * fgrow_winter;
      
    } //End cell overwinter
    
    if(Debug == 1){
      printf("End overwintering. End year. \n");
    }
    
    //Free storage:
    free_vector(oldN,1,Cells);
    free_vector(oldZ,1,Cells); 
    free_vector(Z_begin,1,Cells);
    free_vector(oldB,1,Cells);
    free_vector(oldF,1,Cells);
    free_vector(Percent_Defol,1,Cells);

  }// End year
  
  
  
  //Free storage:
  free_vector(N,1,Cells); 
  free_vector(Z,1,Cells); 
  free_vector(B,1,Cells); 
  free_vector(F,1,Cells); 

  //printf("Succenn:Read in Params...\n");
  // ------------------------------------
  // END ParamNum anningment
  
  
}

