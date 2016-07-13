// --------------------------------Begin ODE system of White model --------------------------------------------//
int fast_odes(double t, const double y[], double dydt[],void *Paramstuff)
{
    //struct STRUCTURE *Params=(struct STRUCTURE *)Paramstuff;
    ParamStruct* ODEParams;
    ODEParams = (ParamStruct*) Paramstuff;
    
    
    double k		= ODEParams->FxdPars[1];
    double V 		= 1/k; //1/(C^2)
    double nu		= ODEParams->FxdPars[2];
    double mu		= ODEParams->FxdPars[3];
    double delta	= ODEParams->FxdPars[6];
    int m			= ((int) ODEParams->FxdPars[7]);  //CK//  the number of exposed classes.  adjustable and fit.  yay!
    double feed     = ODEParams->FxdPars[12];
    double bgrow    = ODEParams->FxdPars[13];
    
    double S0       = ODEParams->FxdPars[14];
    double Fol0     = ODEParams->FxdPars[15]; //Initial Foliage density
    
    
    double N; //total pop size for the foliage model
    
    
    int i;
    
    double hetero = y[m+1]*nu*pow( ((double) (y[0]/S0)), ((double) V));	// save time by doing once (heterogeneity term), where y[m+1] is P, y[0] is S, and V = C^2,
    // ------------------------------------------ ODEs -------------------------------------------- //
    
    
    if(y[0]< 1.0e-12){
        dydt[0]=0;
    }else{
        dydt[0]  = -y[0]*hetero;	 //dS/dt
    }
    
    dydt[1]  = y[0]*hetero - m*delta*y[1]; //dE1/dt
    
    for(i=2; i <= m; i++){  //dEi/dt for i = 2,m
        dydt[i]=m*delta*(y[i-1] -y[i]);
    }
    
    dydt[m+1]  = m*delta*y[m] - mu*y[m+1];  //dP/dt, where y[m] is E_m
    dydt[m+2] = m*delta*y[m];  //dR/dt
    
    //Now do the foliage model:
    for(i=0; i <= m; i++){
        N += y[i]; //calculate the total pop size (S + sum(E))
    }
    
    dydt[m+3] = -feed * N * y[m+3]; // Foliage density
    dydt[m+4] = (y[m+4] * bgrow * (y[m+3] / Fol0)); // Biomass density
	   
    
    return GSL_SUCCESS;
}


// ------------------------------------------  ODE Solver  ----------------------------------------------- //
double ODE_Solver(double t_ode,double t_end,void *Paramstuff,double *y_ode)
{
    int i;
    int status_ode;
    double h_init=1.0e-5;
    
    ParamStruct* ODEParams;
    ODEParams = (ParamStruct*) Paramstuff;
    
    int DIM = ODEParams->FxdPars[7]+5; //The total number of ODEs is the number of E classes, plus 3 because there are also the S, P and R classes + 1 for Foliage + 1 for Biomass
    
    
    const gsl_odeiv2_step_type *solver_ode	= gsl_odeiv2_step_rkf45; // Runge-Kutta Fehlberg (4, 5)
    //const gsl_odeiv_step_type *solver_ode = gsl_odeiv_step_rk4;
    
    // returns pointer to a newly allocated instance of a stepping function of type 'solver_ode' for a system of DIM dimensions //
    gsl_odeiv2_step *step_ode	= gsl_odeiv2_step_alloc(solver_ode, DIM);
    
    //gsl_odeiv_control *tol_ode	= gsl_odeiv_control_standard_new(1.0e-10, 1.0e-5, 1.0, 0.2);
    gsl_odeiv2_control *tol_ode	= gsl_odeiv2_control_standard_new(1.0e-10, 1.0e-5, 1.0, 0.2);
    
    gsl_odeiv2_evolve *evol_ode	= gsl_odeiv2_evolve_alloc(DIM);
    
    //gsl_odeiv_system sys_ode;
    gsl_odeiv2_system sys_ode;
    sys_ode.function  = fast_odes;
    sys_ode.jacobian  = NULL;
    sys_ode.dimension = (size_t)(DIM);
    sys_ode.params	  = ODEParams;
    
    //double y_err[DIM]; double dydt_in[DIM];	double dydt_out[DIM];
    
    // ----------------------------------- Integrate Over Time ------------------------------------ //
    
    
    
    while (t_ode<t_end)	{
        
        
        status_ode = gsl_odeiv2_evolve_apply(evol_ode, tol_ode, step_ode, &sys_ode, &t_ode, t_end, &h_init, y_ode);
        
        //status_ode = gsl_odeiv2_step_apply(step_ode, t_ode, h_init, y_ode, y_err, dydt_in, dydt_out, &sys_ode);
        
        
        
        for (i=0;i<DIM;i++)	{
            
            
            if(y_ode[i]<0) y_ode[i] = 0.0;
            //This next bit holds only because initial cadavers is a fraction in RK45SEIR
            // 		if (y_ode[i]>ODEParams->FxdPars[11])	{  //FxdPars[11] is S0
            // 			//printf("y(%d) TOO LARGE!!\n",i);  
            // 			y_ode[i]=ODEParams->FxdPars[11];
            // 		}
            
        }
        
        //t_ode+=h_init;
        
        
    }
    // -------------------------------------- Clear Memory ----------------------------------------- //
    
    
    gsl_odeiv2_evolve_free(evol_ode);
    gsl_odeiv2_control_free(tol_ode);
    gsl_odeiv2_step_free(step_ode);
    
    return (t_end);
}


