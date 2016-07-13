
float * DistdRK45(double *Pars, float *ModelFI, float *ModelS, double *RandNums, int MaxDays){
    
    FILE *fp;
    ParamStruct ODEParams;
    
    double ratio;
    double z;
    double s;
    
    double zinit;  /* initial pibs, to simulate overwintering */
    double SumE; //sum of the E classes
    int i;
    double t;
    float WM[MaxDays], Suscept[MaxDays];
    double FractI;
    double S0;
    double nubar, mu, k, gamma, delta, feed, bgrow;
    double nubarAdj;
    int MaxStages;
    double S0Adj;
    double TooHigh = 1e10;
    int TooHighOut;
    double tempInt;
    float Z0;
    double Fol0, Fol, Biom0, Biom;
    
    //Store Output:
    static float OUT[3];
    
    int Figures=0;
    if(Figures==1) fp = fopen("FullSEIR.dat","a");
    
    k = Pars[1];
    nubar = Pars[2];
    mu = Pars[3];
    ratio = Pars[5];
    delta = Pars[6];
    
    double StagesTemp = Pars[7] + 0.5;
    int m = ( (int) StagesTemp);   //Number of Stages, rounded off
    
    feed = Pars[8];
    bgrow = Pars[9];
    
    for(i=0;i<MaxDays;i++){
        WM[i] = 0.0;
        Suscept[i] = 0.0;
    }
    
    //Initial density of virus and susceptibles
    S0 = Pars[10]; //Initial population density
    Z0 = Pars[11]; //Initial virus/cadaver density
    s = S0;
    z = Z0;
    
    //Biomass initial density
    Biom0 = Pars[12];
    Biom = Biom0;
    
    //Foliage initial density
    Fol0 = Pars[13];
    Fol = Fol0;
    
    
    //printf("z:%f\t s:%f\t ratio:%f\t FractI:%f\n",z,s,ratio,FractI);
    
    S0Adj = S0;  //Initial density of hosts, saved for later use
    
    int OldTimeIntrvl = -1;  //Useful for checking the time interval action
    
    double y_ode[m+3];
    
    for(i=1;i<=m;i++) { y_ode[i] = 0; } //These are the E classes
    y_ode[0] = s;
    y_ode[m+1] = z;
    y_ode[m+2] = 0; //Dead
    y_ode[m+3] = Fol0; //Foliage
    y_ode[m+4] = Biom0; //Biomass
    
    
    ODEParams.FxdPars[1] = k;
    ODEParams.FxdPars[3] = mu;
    ODEParams.FxdPars[6] = delta;
    ODEParams.FxdPars[7] = (double) m;
    ODEParams.FxdPars[12] = feed;
    ODEParams.FxdPars[13] = bgrow;
    ODEParams.FxdPars[14] = S0;
    ODEParams.FxdPars[15] = Fol0;
    
    double oldt = 0.0;
    int Day = 0;
    
    SumE = 0.0;
    
    double TooLow = 1e-10;
    int Crash = 0;
    //Loop over t
    t = 0; TooHighOut = 0;
    
    //The day upon which Foliage becomes 0 (total defoliation)
    double t_F0 = 0;
    
    while((t<MaxDays)&&(Crash!=1)){
        
        if(Figures==1){
            fprintf(fp,"%f %e %e %e %e %e\n",t,s,z,SumE,Fol,Biom);
            //printf("%f %e %e %e\n",t,s,z,SumE);
        }
        
        nubarAdj = RandNums[Day]; //TEMP
        
        ODEParams.FxdPars[2] = nubarAdj;
        
        ///////////////////  Here is the ODE solver at work ///////////////////
        double t_next = t + 1;
        
        t = ODE_Solver(t,t_next,&ODEParams,y_ode); //Numerical integration from t to t_next, returns (t_end)
        s = y_ode[0];
        z = y_ode[m+1]; //set pop densities equal to ode output
        Fol = y_ode[m+3];
        Biom = y_ode[m+4];
        
        
        // Check for total defoliation
        if(Fol < 1e-9) Fol=0.0;
        
        if(t_F0==0.0){
            if(Fol==0.0) t_F0 = t; //
        }
        
        
        //Sum up exposed classes:
        SumE = 1e-300; //For some horrible reason, this works better than SumE=0;
        for(i=1;i<=m;i++) SumE = SumE + y_ode[i];
        //printf("SumE:%f\n", SumE);
        y_ode[m+2] = 0.0;
        Day++;//update Day number
        
        
        Suscept[Day] = s;
        
        
        if((s<=0.0)&&(SumE<=0.0))
            WM[Day] = 0.0;
        else
            WM[Day] = SumE/(s+SumE);
        
        
        if(WM[Day]>1) WM[Day] = 1.0;
        if(WM[Day]<0) WM[Day] = 0.0;
        
        
    } // t loop 
    
    
    for(i=0;i<Day;i++){
        if((isnan(WM[i])!=1)&&(TooHighOut!=1))
            ModelFI[i] = WM[i];
        else
            ModelFI[i] = 0.0;
        //printf("ModelFI:%f\n",ModelFI[i]);
        
        ModelS[i] = Suscept[i];
        
    }
    
    if(Figures==1){ 
        fclose(fp);
        //exit(1);
    }
    
    //Store output of Foliage and Biomass
    if(Fol < 1e-10) Fol = 0;
    if(Biom < 1e-10) Biom = 0;
    
    OUT[0] = Fol;
    OUT[1] = t_F0; //this will be 0.0 if Fol > 0
    OUT[2] = Biom;

    
    return OUT;
    
}  //The End






