

#define IntMax 600000
// tau divided by the step size


int PrintMe = 0; //whether to print output to a file


double F(double sus,double nubar, double path, double sus0, double k)
 {
    return( - sus * nubar * path * pow(sus/sus0,1/k));
  }

double G(double sus,double nubar, double path, double sus0, double k, double inf, double gamma)
  {
   return(sus * nubar * path * pow(sus/sus0,1/k) - gamma*inf);
  }

double G2(double inf1, double gamma, double inf2)
{
    return(gamma*inf1 - gamma*inf2);
}



double H(double inf, double gamma, double path, double mu)
  {
   return(gamma*inf - mu*path);
  }


double Distd1(double *Pars,float *Covars){


 	//This version assumes that Covars[2][1] is the actual virus pop, NOT the fractinf 
 double initz;
 int flag = 1;
 int Plt;
 //int MaxPlts = 5;

 double ratio;
 int Split;
 double Hch;
 double maxData;
 int plt;
 FILE *fp,*fp2,*fp3;
 double maxt;    // number of days to run for
 //double l,l0,l1,l2,l3;
 double n,n0,n1,n2,n3; //Dead
 //double p,p0,p1,p2,p3;
 double q,q0,q1,q2,q3;
 double h;  /* time step */

 double z;
 int ind; //indicates where you are in storage arrays
 double *Intx;
 double *SStor, *PStor, *nuStor;
 double Stmp, Ptmp, nutmp;
 double l, m[100], p;
 double l0, l1, l2, l3;
 double m0[100], m1[100], m2[100], m3[100];
 double p0, p1, p2, p3;
 //double l0tmp, l1tmp, l2tmp, l3tmp;
 //double p0tmp, p1tmp, p2tmp, p3tmp;
 
 double s, inf[100];
 /* *initx; */
  /* z is pibs, and s is suscept.'s */
 /* Intx is sum over the suscepts. times nu times pibs z, for storage */
 double zinit;  /* initial pibs, to simulate overwintering */
 double Dead,y,oldy;      /* infect.'s - no need to store'em */
 int i;
 double t;
 int tau, tau0;
 float WM[1000], Suscept[1000];
 //double *WM, *Suscept;
 int week;
 double FractI;
 double S0;
 double nubar, mu, k, gamma;
 double nubarAdj; 
 int MaxStages;
 int CorrlnIntrvl, TimeIntrvl;
 double DataDay;
 double S0Adj;
 double TooLow = -1;
 double TooHigh = 1e10;
 int TooHighOut;
 double tempInt;
 int OutputCount;
 int ReDo, ReDoCount = 0;
 int Bang, FirstCrash;
 float Z0;
 double SExact, IExact, PExact, dummy1, dummy2, dummy3;

 //Intx = dvector(0,IntMax-1);  //allocates memory, more than needed
 SStor = dvector(0,IntMax-1); 
 PStor = dvector(0,IntMax-1); 
 nuStor = dvector(0,IntMax-1);


  k = Pars[2];  //model parameter: Squared CV of dist'n of transmission rates
   nubar = Pars[1];   //average transmission rate


  mu = Pars[3];

  h = Pars[7];   //DDE time step

  CorrlnIntrvl = Pars[8];  //Correlation Interval

  MaxStages = Pars[9];   //Number of Realizations

  ratio = pow(10,Pars[5]); //size of early stage larva to size of later stage larva, in # of virus particles



ReDo = 1; h = Pars[7];
TooHighOut = 0;
Plt=Pars[0]; 


	for(i=0;i<20;i++){
		WM[i] = 0.0;
		Suscept[i] = 0.0;
	}

	//tau0 = Pars[6];

	//tau = Pars[6]/h;  //time delay is tau, divide by h to get how much memory you need
	gamma = Pars[6];

    S0 = Covars[1]; //Initial population size
 	FractI = Covars[2]; //Initial fraction infected
	Split = 0;  //Whether hatch time is spread out, yes or no


if(FractI<0){
 FractI = pow(10,Pars[11+Plt]);
}


	 maxt = 56.0;  

   if(Split)
    Hch = 0.5;  //If there is some spread, half hatches in first week
   else
    Hch = 1.0;  //No spread, everyone hatches at once.

   //z = S0*FractI*Hch*ratio;  //Initial density of virus
   
   z = Covars[2]; //THIS IS THE ONLY (MAIN?) DIFF OF DDE_SSE9B from DDE_SSE9
   //printf("S0:%f z:%f nubar:%f mu:%f \n",S0,z,nubar,mu);
   Z0 = z;
   //s = S0*(1-FractI);  //Initial density of hosts
   s = S0; //this is for graphs for Emma's paper
   for(i=0;i<MaxStages;i++)
	inf[i] = 0;

   S0Adj = s;  //Initial density of hosts, saved for later use
   initz = z;  //Saving initial virus density, not sure why?
   //printf("S0Adj:%f z:%f k:%f\n",S0Adj,z,k);
   //getc(stdin);

 

  week = 1;  //Initializing week no, need it for data storage

  flag = 1;  //Says whether or not additional hatch has occurred

  y= 0.0;  //Initializing number dead per unit time
  Dead = 0.0;  //Initializing cumulative number dead per week
  int OldTimeIntrvl = -1;  //Useful for checking the time interval action
  DataDay = 0.0; 

  OutputCount = 1; 

double Oldnubar = -1;
Bang = 0; FirstCrash = 1;
double StageAdjust;

if(SuperVerbose==1){
	fp = fopen("Test.dat","w");
	fclose(fp);
}

//printf("s:%f z:%f \n",s,z);

ind = 0;


/*
	printf("%f %f %f ",t,s,z);
	for(i=0;i<MaxStages;i++)
		printf("%f ",inf[i]);
	printf("\n");
*/

for(t=h;t<=maxt;t+=h){  //time loop

  /*
	printf("%f %f %f ",t,s,z);
	for(i=0;i<MaxStages;i++)
		printf("%f ",inf[i]);
	printf("\n");
  */
	
		TimeIntrvl = (t/CorrlnIntrvl);  //TimeIntrvl indexes the Random Numbers - you want it to truncate the floats

	
		//nubarAdj = (RandNums[1][Plt][TimeIntrvl]);
		//nutmp = nuStor[ind];
		nubarAdj = nubar;
		nutmp = nubar;


    //DDE solver: method of steps using 4th order R-K
		
	//l0tmp = l0[ind]; l1tmp = l1[ind]; l2tmp = l2[ind]; l3tmp = l3[ind]; 
	//p0tmp = p0[ind]; p1tmp = p1[ind]; p2tmp = p2[ind]; p3tmp = p3[ind];
	//Stmp = SStor[ind]; Ptmp = PStor[ind];

    l0 = F(s,nubarAdj,z,S0Adj,k);
    
	m0[0] = G(s,nubarAdj,z,S0Adj,k,gamma,inf[0]);
	for(i=1;i<MaxStages;i++)
		m0[i] = G2(inf[i-1],gamma,inf[i]);
    p0 = H(gamma,inf[MaxStages-1],mu,z);
		
    l1 = F(s+0.5*l0*h,nubarAdj,z+0.5*p0*h,S0Adj,k);
	m1[0] = G(s+0.5*l0*h,nubarAdj,z+0.5*p0*h,S0Adj,k,inf[0]+0.5*m0[0]*h,gamma);
	for(i=1;i<MaxStages;i++)
		m1[i] = G2(inf[i-1]+0.5*m0[i-1]*h,gamma,inf[i]+0.5*m0[i]*h);
	p1 = H(inf[MaxStages-1]+0.5*m0[MaxStages-1]*h,gamma,z+0.5*p0*h,mu);
	
    l2 = F(s+0.5*l1*h,nubarAdj,z+0.5*p1*h,S0Adj,k);
	m2[0] = G(s+0.5*l1*h,nubarAdj,z+0.5*p1*h,S0Adj,k,inf[0]+0.5*m1[0]*h,gamma);
	for(i=1;i<MaxStages;i++)
		m2[i] = G2(inf[i-1]+0.5*m1[i-1]*h,gamma,inf[i]+0.5*m1[i]*h);
    p2 = H(inf[MaxStages-1]+0.5*m1[MaxStages-1]*h,gamma,z+0.5*p1*h,mu);
	
    l3 = F(s+l2*h,nubarAdj,z+p2*h,S0Adj,k);
	m3[0] = G(s+l2*h,nubarAdj,z+p2*h,S0Adj,k,inf[0]+m2[0]*h,gamma);
	for(i=1;i<MaxStages;i++)
		m3[i] = G2(inf[i-1]+m2[i-1]*h,gamma,inf[i]+m2[i]*h);
	p3 = H(inf[MaxStages-1]+m2[MaxStages-1]*h,gamma,z+p2*h,mu);
   
    
    l = (l0 + 2*l1 + 2*l2 + l3)/6.0;
	for(i=0;i<MaxStages;i++)
		m[i] = (m0[i] + 2*m1[i] + 2*m2[i] + m3[i])/6.0;
    p = (p0 + 2*p1 + 2*p2 + p3)/6.0;



	int Stop = 0;
	/*
        if((fabs(l0)>TooHigh)) Stop = 1;
	   if((fabs(l1)>TooHigh)) Stop = 1;
	      if((fabs(l2)>TooHigh)) Stop = 1;  
		 if((fabs(l3)>TooHigh)) Stop = 1; 
		    if(Stop!=1)
		      s += l*h;
		    else
		      s = 0;
	*/

           s += l*h;
	

	if(s<1e-50) s = 0.0;
	z += p*h;
	for(i=0;i<MaxStages;i++)
	     inf[i] += m[i]*h;
	

if(SuperVerbose==1){
	fp = fopen("Test.dat","a");
	fprintf(fp,"%f %f %f ",t,s,z);
	for(i=0;i<MaxStages;i++)
		fprintf(fp,"%f ",inf[i]);
	fprintf(fp,"\n");
	//printf("%f %f %f %f\n",t,s,inf[MaxStages-1],z);
	//printf("%f %f %f %f %f\n",t,s,y,z,Dead);
	//fprintf(fp,"%f %f %f\n",t,y/(s+y),s);
	fclose(fp);
	//getc(stdin); 
}



   


 }    // t loop 




 //free_dvector(Intx,0,IntMax-1);  //freeing up the memory

 free_dvector(SStor,0,IntMax-1);
 free_dvector(PStor,0,IntMax-1);
 free_dvector(nuStor,0,IntMax-1);


 //if(SuperVerbose==1){ exit(1);}


    double answer = 1.0 - (s/S0Adj);
    //printf("answer:%f s:%f S0Adj:%f\n",answer,s,S0Adj);
	if((answer<0)||(answer>1)) answer = -1;
	//return 1.0 -(s/S0Adj);
	return answer;


}  //The End






