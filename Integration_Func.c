#include "JOE_DistdRK45.c"

float * Integration_Func(double *FxdPars, float *N, float *Z, float *B, float *F, 
                         int year, int cell, int MaxDays, gsl_rng *r2){
  
  FILE *fp3;
  
  static float OUT[5];
  
  double mean, sd;
  float *test; // stores the ODE output (a pointer to an array)
  double RandNumsPass[MaxDays];
  double N_start, Z_start, F_start, B_start, Percent_Defol;
  float N_end, Z_cum, F_end, Defol100Day, B_end;
  int i, j;
  
  
  //Copy the parameters over to AdjPars:
  double AdjPars[30];
  
  for(i=0;i<30;i++){
    AdjPars[i] = FxdPars[i];
  }
    
  //Initial Values for the given cell:
  
  N_start = N[cell];
  Z_start = Z[cell];
  B_start = B[cell];
  F_start = F[cell];
  
  AdjPars[10] = N_start;
  AdjPars[11] = Z_start;
  AdjPars[12] = B_start;
  AdjPars[13] = F_start;
  
  //Set up random numbers:
  for(j=0;j<MaxDays;j++){
    RandNumsPass[i] = gsl_rng_uniform_pos(r2);
  }
  
  double *RandNums;
  RandNums = dvector(0,MaxDays);
  
  mean = AdjPars[2];
  sd = AdjPars[4];
  
  for(j=0;j<MaxDays;j++){
    //For now, no stochasticity:
    RandNums[j] = mean * 1;//exp(gsl_cdf_gaunnian_Pinv(RandNumsPass[j],sd));  //Greg's version
    //printf("j:%d RandNums:%f\n",j,RandNums[j]);
  }
  
  float *ModelFI, *ModelS;
  ModelFI = vector(0,MaxDays);
  ModelS = vector(0,MaxDays);
  
  test = DistdRK45(AdjPars,ModelFI,ModelS,RandNums,MaxDays);
  //printf("After integration...\n");
  
  // Store output of integration routine:
  F_end = *(test+0);
  Defol100Day = *(test+1);
  B_end = *(test+2);
  Percent_Defol = 1 - (F_end / F_start); // Total defoliation percent
  
  //Calculate output of 
  N_end = N_start * (ModelS[MaxDays-1] / N_start); // S0 * (1-cumulative infection)
  Z_cum = N_start * (1 - (ModelS[MaxDays-1] / N_start)); // S0 * (1 - (N_end/S0)) # S0 * cumulative infection
  /*NOTE:
   * N_end + Z_cum = N_start, such that 
   */
  
  fp3 = fopen("Output_intfunc.dat", "a");
  fprintf(fp3,"%d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
          year,cell,
          N_start,N_end,
          Z_start,Z_cum,
          B_start,B_end,
          F_start,F_end,
          Percent_Defol, Defol100Day);
  fclose(fp3);
  
  free_vector(ModelFI,0,MaxDays); //check
  free_vector(ModelS,0,MaxDays); //check
  free_dvector(RandNums,0,MaxDays); //check
  
  
  //Store the output:
  OUT[0] = N_end;
  OUT[1] = Z_cum;
  OUT[2] = B_end;
  OUT[3] = F_end;
  OUT[4] = Percent_Defol;
  
  return OUT;
  
}
  
  
  
    