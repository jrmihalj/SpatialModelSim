
float * Integration_Func(double *FxdPars, float *N, float *Z, float *B, float *F, int year, int cell, int MaxDays, gsl_rng *r2){
  
  FILE *fp2;
  
  static float OUT[4];
  
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
  fp2 = fopen("Output.dat", "a");
  
  fprintf(fp2,"%d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
          year,cell,
          N_start,N_end,
          Z_start,Z_cum,
          B_start,B_end,
          F_start,F_end,
          Percent_Defol, Defol100Day);
  
  fclose(fp2);
  
  free_vector(ModelFI,0,MaxDays); //check
  free_vector(ModelS,0,MaxDays); //check
  free_dvector(RandNums,0,MaxDays); //check
  
  
  //Store the output:
  OUT[0] = N_end;
  OUT[1] = Z_cum;
  OUT[2] = B_end;
  OUT[3] = F_end;
  
  return OUT;
  
}
  
  
  
  
  
  /* UNUSED:
   
   
   
   for(zz = Lower_Z; zz <= Upper_Z; zz += Jump_Z){
   
   if(Fixed == 0){
   AdjPars[10] = zz;
   }else{
   AdjPars[10] = z_fixed; //really z initial
   zz = Upper_Z;
   }
   
   
   for(nn = Lower_N; nn <= Upper_N; nn += Jump_N){
   
   if(Fixed==0){
   AdjPars[9] = exp(nn);
   }else{
   AdjPars[9] = n_fixed;
   nn = Upper_N;
   }
   
   N_start = AdjPars[9];
   Z_start = AdjPars[10];
   F_start = AdjPars[11];
   B_start = AdjPars[13];
   
   first = 1;
   
   for(i=1;i<=reps;i++){
   
   for(j=0;j<MaxDays;j++){
   RandNumsPass[i] = gsl_rng_uniform_pos(r2);
   }
   
   double *RandNums;
   RandNums = dvector(0,MaxDays);
   
   mean = AdjPars[2];
   sd = AdjPars[4];
   
   for(j=0;j<MaxDays;j++){
   RandNums[j] = mean * 1;//exp(gsl_cdf_gaunnian_Pinv(RandNumsPass[j],sd));  //Greg's version
   //printf("j:%d RandNums:%f\n",j,RandNums[j]);
   }
   
   float *ModelFI, *ModelS;
   ModelFI = vector(0,MaxDays);
   ModelS = vector(0,MaxDays);
   
   test = DistdRK45(AdjPars,ModelFI,ModelS,RandNums,MaxDays,yearly,first,Fixed);
   //printf("After integration...\n");
   
   Fol_end = *(test+0);
   Defol100Day = *(test+1);
   Biom_end = *(test+2);
   Percent_Defol = 1 - (Fol_end / F_start); // Total defoliation percent
   
   if(Fixed==0){
   fp = fopen("Output.dat", "a");
   }else{
   fp = fopen("Output_eqm.dat", "a");
   }
   
   //year, rep, FI_end, N_end, F_start, Fol_end, 100DefolDay, B_start, Biom_end, Z_n, N_n, Z0, N0
   fprintf(fp,"%d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
   1,i,ModelFI[MaxDays-1],ModelS[MaxDays-1],F_start,Fol_end,Defol100Day,B_start,Biom_end,Z_start,N_start,zz,exp(nn));
   
   fclose(fp);
   
   if(yearly){
   N_end = N_start * (ModelS[MaxDays-1] / N_start); // S0 * (1-cumulative infection)
   Z_cum = N_start * (1 - (ModelS[MaxDays-1] / N_start)); // S0 * (1 - (N_end/S0)) # S0 * cumulative infection
   
   }
   
   free_vector(ModelFI,0,MaxDays); //check
   free_vector(ModelS,0,MaxDays); //check
   free_dvector(RandNums,0,MaxDays); //check
   
   first=0;
   
   /* Now do the longitudinal simulation */
   
   /* 
    * 
    * 
    * STILL UNUSED:
    * 
    if(yearly==1){
    
    int remain;
    int y;
    for(y=2;y<=MaxYear;y++){
    
    if(N_end > 0.0){
    
    float *ModelFI, *ModelS;
    ModelFI = vector(0,MaxDays);
    ModelS = vector(0,MaxDays);
    
    //Set up RandNums
    double *RandNums;
    RandNums = dvector(0,MaxDays);
    
    mean = AdjPars[2];
    sd = AdjPars[4];
    
    for(j=0;j<MaxDays;j++){
    RandNums[j] = mean * 1;//exp(gsl_cdf_gaunnian_Pinv(RandNumsPann[j],sd));  //Greg's version
    }
    //End RandNums set up
    
    //Calculate starting conditions for N w/ predator model:
    AdjPars[9] = N_end * reprod; //* (1 - ((a_pred*b_pred*N_start)/(pow(b_pred,2)+pow(N_start,2))));
    N_start = AdjPars[9];
    //Calculate starting conditions for Z:
    AdjPars[10] = Z_cum * winter + gamma * Z_start; //This Z_start is from the previous generation, so Z_n
    Z_start = AdjPars[10]; //Now Z_start has been over-written, so Z_n+1 (i.e. the current generation)
    
    //Calculate starting conditions for Biomass:
    AdjPars[13] = bgrow_winter * Biom_end * (1 - (a_biom * Percent_Defol / (c_biom + Percent_Defol)));
    B_start = AdjPars[13];
    
    //Calculate starting conditions for Foliage:
    AdjPars[11] = Fol_end + Biom_end * fgrow_winter;
    F_start = AdjPars[11];
    
    
    //Run integration:
    test = DistdRK45(AdjPars,ModelFI,ModelS,RandNums,MaxDays,yearly,first,Fixed);
    
    Fol_end = *(test+0);
    Defol100Day = *(test+1);
    Biom_end = *(test+2);
    Percent_Defol = 1 - (Fol_end / F_start); // Total defoliation percent
    
    if(Fixed==0){
    fp = fopen("Output.dat", "a");
    }else{
    fp = fopen("Output_eqm.dat", "a");
    }
    
    //year, rep, FI_end, N_end, Fol_end, 100DefolDay, Biom_end, Z_n, N_n, Z0, N0
    fprintf(fp,"%d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
    y,i,ModelFI[MaxDays-1],ModelS[MaxDays-1],F_start,Fol_end,Defol100Day,B_start,Biom_end,Z_start,N_start,zz,exp(nn));
    
    fclose(fp);
    
    N_end = N_start * (ModelS[MaxDays-1] / N_start); // S0 * cumulative infection
    Z_cum = N_start * (1 - (ModelS[MaxDays-1] / N_start)); // S0 * (1 - (N_end/S0)) # S0 * cumulative infection
    
    free_vector(ModelFI,0,MaxDays); //check
    free_vector(ModelS,0,MaxDays); //check
    free_dvector(RandNums,0,MaxDays); //check
    
    }
    /*
     else{ //if N_end <=0.0
     if(Fixed==0){
     fp = fopen("Output.dat", "a");
     }else{
     fp = fopen("Output_eqm.dat", "a");
     }
     for(remain=y;remain<=MaxYear;remain++){
     fprintf(fp, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",remain,i,0.0,0.0,0.0,0.0,zz,exp(nn));
     }
     y=MaxYear;
     fclose(fp);
     } */
    
    /* STILL UNUSED:
     
     
    }//for each year
     
     
    }//if yearly==1
     
     
     
   }//for each rep (stochastic sims)
     
   }//for each S0
     
   }//for each FI0
     
     
     
     
     return 0;
     
     }
     
     End unused*/ 
    
    