

#include "time.h"
#include "string.h"
#include "stdlib.h"

double PRD = 1e4;  //This is the "periodicity" of stable sub-populations.  The fraction of sub-pops that are stable = 1/PRD
//double PROB = 0.25; //Probability of being a non-inducible forest type
double PI = 3.14159265;

//#include "nr.h"
#include "nrutil.h" 
#include "nrutilGD2.c"

#include "ran2.c"
//#include "gammln.c"
//#include "poidev.c"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>

gsl_rng *r2;


#include "Zbrentd2.c"
#include "FX.C"

//#include "ran1.c"
//#include "gasdev.c"



#include "pearsnD.c"

#include "CorrldRand.c"

#include "time.h"

int SuperVerbose = 0;

#include "Distd2.c"


//long seed = -5;

void Initialize(double **N, double **Z, double **D, double **ForestType, double InitN, double InitZ, double InitD, int MaxRow, int MaxCol,float IndProb){

	int i,j;
	long seed;

	int StartPoint = ceil((double) MaxRow/2);
	
	//seed = -1;
	
	//srand((unsigned)time(NULL)); // initialize ran # generator
 	//long idum= -rand();
	//seed = idum;

	//double rnum = ran2(&seed);
        double rnum;

	for(i=1;i<=MaxRow;i++)
		for(j=1;j<=MaxCol;j++){

		  //N[i][j] = Z[i][j] = 0; D[i][j] = 0.01; //really this should be 0, so that it can, if it wants, find its way to some value above InitD;
		  N[i][j] = InitN; //Here we assume that the gypsy moth invasion has already occurred
		  Z[i][j] = InitZ;
		  D[i][j] = 0.01;  //But the induced defenses start at a very low level.
		  rnum = gsl_rng_uniform(r2);
		  //rnum = ran2(&seed);


		  if(rnum<IndProb){
				ForestType[i][j] = 0;
		  }	else
				ForestType[i][j] = 1;
		} //for j
		//The next 3 lines are kind of pointless, but they allow for the possibility that you want gypsy moths and virus at only 3 initial points
		N[1][StartPoint] = InitN; Z[1][StartPoint] = InitZ;
		N[1][StartPoint+1] = InitN; Z[1][StartPoint+1] = InitZ;
		N[1][StartPoint-1] = InitN; Z[1][StartPoint-1] = InitZ;
} //end of Initialize


//%CALCULATE NORMALIZING CONSTANT FOR BALLOONING%
double InitBalloon(double kernelCOEF, int MaxRow, int MaxCol, float RowScale, float ColScale){

	int x,y;
	float Dist, DispFract;


	//calculate normalizing constant
	double truncate = 2*RowScale; //max dispersal allowed = 2 grid cells (prob of even ballooning 1 km is almost 0)
	int middleRow = ceil(truncate/RowScale)+1;
	int middleCol = ceil(truncate/RowScale)+1;

	double SumBalloonFract = 0.0;

	for(x=1;x<=2*ceil(truncate/RowScale)+1;x++)
		for(y=1;y<=2*ceil(truncate/ColScale)+1;y++){

			Dist = sqrt(  pow(((double)(x-middleRow)*RowScale),2) + pow((double) ((y-middleCol)*ColScale),2)  );
			if(Dist<=truncate){
				DispFract = (2/(kernelCOEF*sqrt(2*3.1416)))*exp(-(pow((Dist/(2*kernelCOEF/1000)),2)));
			}else{
				DispFract = 0;}
			SumBalloonFract += DispFract;
			
		}
	//end norm const calculation
	return SumBalloonFract;

} //InitBalloon


//%CALCULATE NORMALIZING CONSTANT FOR LONG-DISTANCE DISPERSAL%
double InitLDD(double LongDistParam, double LongDistFract, int MaxRow, int MaxCol, float RowScale, float ColScale, double truncateProb){

	int x,y;
	float Dist, DispFract;
	
////////////////////////////// EXPONENTIAL KERNEL ///////////////////////////
	//calc normalizing const
	float truncateLDD = -LongDistParam*log(truncateProb*LongDistParam); //max km for LDD
	int middleRow = ceil(truncateLDD/RowScale)+1;
	int middleCol = ceil(truncateLDD/ColScale)+1;

	double SumLDDFract = 0.0;

	for(x=1;x<=2*ceil(truncateLDD/RowScale)+1;x++){
		for(y=1;y<=2*ceil(truncateLDD/ColScale)+1;y++){

			Dist = sqrt(  pow((double)((x-middleRow)*RowScale),2) + pow((double)((y-middleCol)*ColScale),2)  );
			if(Dist<=truncateLDD){
				DispFract = (1/LongDistParam)*exp(-(Dist/LongDistParam));
			}else{
				DispFract = 0;}
			SumLDDFract += DispFract;
			
	
		} 
        }
		
		
////////////////////////////// NORMAL KERNEL ///////////////////////////
	//calc normalizing const
/*	float truncateLDD = 2*LongDistParam*sqrt(log(2/(truncateProb*LongDistParam*sqrt(2*3.1416)))); 
	int middleRow = ceil(truncateLDD/RowScale)+1;
	int middleCol = ceil(truncateLDD/ColScale)+1;

	double SumLDDFract = 0.0;

	for(x=1;x<=2*ceil(truncateLDD/RowScale)+1;x++)
		for(y=1;y<=2*ceil(truncateLDD/ColScale)+1;y++){

			Dist = sqrt(  pow((double (x-middleRow)*RowScale),2) + pow((double (y-middleCol)*ColScale),2)  );
			if(Dist<=truncateLDD){
				DispFract = (2/(LongDistParam*sqrt(2*3.1416)))*exp(- pow((Dist/(2*LongDistParam)),2));
			}else{
				DispFract = 0;}
			SumLDDFract += DispFract;
			
 		} */

/////////////////////////////////////////////////////////////////////////

	//end norm const calc
	return SumLDDFract;

} //InitLDD

void DispersalFuncOttar(double **N, double **Z, double SumBalloonFract, double kernelCOEF, double DispFract, int MaxRow, int MaxCol, int RowScale, int ColScale){

	double **NIncoming, **NOutgoing;
	double **ZIncoming, **ZOutgoing;
	int i,j,k,l;
	float Dist;
	double MovingN, MovingZ;

	NIncoming = dmatrix(1,MaxRow,1,MaxCol);
	NOutgoing = dmatrix(1,MaxRow,1,MaxCol);
	ZIncoming = dmatrix(1,MaxRow,1,MaxCol);
	ZOutgoing = dmatrix(1,MaxRow,1,MaxCol);

	double truncate = 2*RowScale;

		for(k=1;k<=MaxRow;k++)
			for(l=1;l<=MaxCol;l++){
				NOutgoing[k][l] = 0;
				ZOutgoing[k][l] = 0;
			}

		MovingN = DispFract*N[1][1];
		NIncoming[1][2] = NOutgoing[1][1] = MovingN;
		MovingN = DispFract*N[1][2];
		NIncoming[1][1] = NOutgoing[1][2] = MovingN;

		N[1][1] += NIncoming[1][1] - NOutgoing[1][1];
		N[1][2] += NIncoming[1][2] - NOutgoing[1][2];

		MovingZ = DispFract*Z[1][1];
		ZIncoming[1][2] = ZOutgoing[1][1] = MovingZ;
		MovingZ = DispFract*Z[1][2];
		ZIncoming[1][1] = ZOutgoing[1][2] = MovingZ;

		Z[1][1] += ZIncoming[1][1] - ZOutgoing[1][1];
		Z[1][2] += ZIncoming[1][2] - ZOutgoing[1][2];

/*
		for(i=1;i<=MaxRow;i++)
			for(j=1;j<=MaxCol;j++){
				NIncoming[i][j] = 0;
				ZIncoming[i][j] = 0;
				for(k=1;k<=MaxRow;k++)
					for(l=1;l<=MaxCol;l++){

							MovingN = DispFract*N[k][l];
							NIncoming[i][j] += MovingN;
							NOutgoing[k][l] += MovingN;

							MovingZ = DispFract*Z[k][l];
							ZIncoming[i][j] += MovingZ;
							ZOutgoing[k][l] += MovingZ;
								
					}
			}

		for(i=1;i<=MaxRow;i++)
			for(j=1;j<=MaxCol;j++){
				
				N[i][j] += NIncoming[i][j] - NOutgoing[i][j];
				if(N[i][j]<0)
					N[i][j] = 0;

				//Z[i][j] += ZIncoming[i][j] - ZOutgoing[i][j];
				if(Z[i][j]<0)
					Z[i][j] = 0;
				
			}
*/


	free_dmatrix(NIncoming,1,MaxRow,1,MaxCol);
	free_dmatrix(NOutgoing,1,MaxRow,1,MaxCol);
	free_dmatrix(ZIncoming,1,MaxRow,1,MaxCol);
	free_dmatrix(ZOutgoing,1,MaxRow,1,MaxCol);

} //DispersalFunc


//%CALCULATE FRACTION TO MOVE INTO/OUT OF EACH POPULATION BY BALLOONING%
void DispersalFunc(double **N, double **Z, double SumBalloonFract, double kernelCOEF, double LongDistFract, int MaxRow, int MaxCol, float RowScale, float ColScale){

	double **NIncoming, **NOutgoing;
	double **ZIncoming, **ZOutgoing;
	int i,j,k,l;
	float Dist, DispFract;
	double MovingN, MovingZ;

	NIncoming = dmatrix(1,MaxRow,1,MaxCol);
	NOutgoing = dmatrix(1,MaxRow,1,MaxCol);
	ZIncoming = dmatrix(1,MaxRow,1,MaxCol);
	ZOutgoing = dmatrix(1,MaxRow,1,MaxCol);

	double truncate = 2*RowScale;

		for(k=1;k<=MaxRow;k++)
			for(l=1;l<=MaxCol;l++){
				NOutgoing[k][l] = 0;
				ZOutgoing[k][l] = 0;
			}

		for(i=1;i<=MaxRow;i++)
			for(j=1;j<=MaxCol;j++){
				NIncoming[i][j] = 0;
				ZIncoming[i][j] = 0;
				for(k=1;k<=MaxRow;k++)
					for(l=1;l<=MaxCol;l++){

							Dist = sqrt(  pow(((double)(i-k)*RowScale),2) + pow(((double) (j-l)*ColScale),2)  );
							if(Dist<=truncate){
								DispFract = (2/(kernelCOEF*sqrt(2*3.1416)))*exp(-(pow((Dist/(2*kernelCOEF/1000)),2)));
							}else{
								DispFract = 0;}

							MovingN = (DispFract/SumBalloonFract)*(1-LongDistFract)*N[k][l];
							NIncoming[i][j] += MovingN;
							NOutgoing[k][l] += MovingN;

							MovingZ = (DispFract/SumBalloonFract)*(1-LongDistFract)*Z[k][l];
							ZIncoming[i][j] += MovingZ;
							ZOutgoing[k][l] += MovingZ;
								
					}
			}

		for(i=1;i<=MaxRow;i++)
			for(j=1;j<=MaxCol;j++){
				
				N[i][j] += NIncoming[i][j] - NOutgoing[i][j];
				if(N[i][j]<0)
					N[i][j] = 0;

				Z[i][j] += ZIncoming[i][j] - ZOutgoing[i][j];
				if(Z[i][j]<0)
					Z[i][j] = 0;
				
			}


	free_dmatrix(NIncoming,1,MaxRow,1,MaxCol);
	free_dmatrix(NOutgoing,1,MaxRow,1,MaxCol);
	free_dmatrix(ZIncoming,1,MaxRow,1,MaxCol);
	free_dmatrix(ZOutgoing,1,MaxRow,1,MaxCol);

} //DispersalFunc

//%CALCULATE FRACTION TO MOVE INTO/OUT OF EACH POPULATION BY LONG-DISTANCE DISPERSALS%
void DispersalFunc2(double **N, double **Z, double SumLDDFract, double LongDistParam, double LongDistFract, int MaxRow, int MaxCol, int RowScale, int ColScale, double truncateProb){

	double **NIncoming, **NOutgoing;
	double **ZIncoming, **ZOutgoing;
	int i,j,k,l;
	float Dist, DispFract;
	double MovingN, MovingZ;
	
	NIncoming = dmatrix(1,MaxRow,1,MaxCol);
	NOutgoing = dmatrix(1,MaxRow,1,MaxCol);
	ZIncoming = dmatrix(1,MaxRow,1,MaxCol);
	ZOutgoing = dmatrix(1,MaxRow,1,MaxCol);

	float truncateLDD = -LongDistParam*log(truncateProb*LongDistParam); //for Exp kernel //max dispersal is the dist with probability of 10^(-6)
	//float truncateLDD = 2*LongDistParam*sqrt(log(2/(truncateProb*LongDistParam*sqrt(2*3.1416)))); //for Norm kernel
	
		for(k=1;k<=MaxRow;k++)
			for(l=1;l<=MaxCol;l++){
				NOutgoing[k][l] = 0;
				ZOutgoing[k][l] = 0;
			}

		for(i=1;i<=MaxRow;i++)
			for(j=1;j<=MaxCol;j++){
				NIncoming[i][j] = 0;
				ZIncoming[i][j] = 0;
				for(k=1;k<=MaxRow;k++)
					for(l=1;l<=MaxCol;l++){

							Dist = sqrt(  pow((double)((i-k)*RowScale),2) + pow((double)((j-l)*ColScale),2)  );
							if(Dist<=truncateLDD){
								DispFract = (1/LongDistParam)*exp(-(Dist/LongDistParam)); // Exponential kernel
								//DispFract = (2/(LongDistParam*sqrt(2*3.1416)))*exp(- pow((Dist/(2*LongDistParam)),2)); // Normal kernel
							}else{
								DispFract = 0;}

							MovingN = (DispFract/SumLDDFract)*LongDistFract*N[k][l];
							NIncoming[i][j] += MovingN;
							NOutgoing[k][l] += MovingN;

							MovingZ = (DispFract/SumLDDFract)*LongDistFract*Z[k][l];
							ZIncoming[i][j] += MovingZ;
							ZOutgoing[k][l] += MovingZ;
					}
			}

		for(i=1;i<=MaxRow;i++)
			for(j=1;j<=MaxCol;j++){
				
				N[i][j] += NIncoming[i][j] - NOutgoing[i][j];
				if(N[i][j]<0)
					N[i][j] = 0;

				Z[i][j] += ZIncoming[i][j] - ZOutgoing[i][j];
				if(Z[i][j]<0)
					Z[i][j] = 0;

				
			}


	free_dmatrix(NIncoming,1,MaxRow,1,MaxCol);
	free_dmatrix(NOutgoing,1,MaxRow,1,MaxCol);
	free_dmatrix(ZIncoming,1,MaxRow,1,MaxCol);
	free_dmatrix(ZOutgoing,1,MaxRow,1,MaxCol);

} //DispersalFunc2

//%POPULATION DYNAMICS%
void EpiDemog(double **N, double **Z, double **D, double **ForestType, int MaxRow, int MaxCol, float *DemogParams, float *EpiParams, double *GaussianNoise){

	int i,j;
	int Pop;
	double StateVars[2];
	float predation;
	double Survivors;
	double FractInf, N2, Z2, D2;
	double cx = 0.0, dx = 1e10, tol = 1e-60;

	float lambda = DemogParams[0];
	float phi = DemogParams[1];
	float gamma = DemogParams[2];
	float a = DemogParams[3];
	float b = DemogParams[4];
	float alpha = DemogParams[5];
	float theta = DemogParams[6];
	float D0 = DemogParams[7];
	float beta = DemogParams[8]; 
        double sigma = DemogParams[9];

	double k;
        int Stop = 0;
	double junk;
	double Pars[10]; //Plt, k, nubar, mu, ?, log10(ratio), tau, h, CorrlnIntrvl, MaxTrials
	float Covars[10];

	Pars[0] = 1.0; //Plt
	Pars[1] = 1.0; //nubar - scaled away
	double nubar = Pars[1]; 
	Pars[2] = 0.015; //baseline k - this is changed down below
	Pars[3] = 0.39; //Emma's value
	Pars[4] =  -1; //sigma
	Pars[5] = log10(2e-3); //ratio, but should this be 1?
        Pars[5] = 1;
	
	Pars[7] = 0.01; //h
	Pars[8] = 1;  //CorrlnIntrvl
	Pars[9] = 25.0; //25; //15; //number of stages, let's say
	Pars[6] = Pars[9]/12.0; //56.0/10.0; //25.0/12.0; //12; //Emma tau
	double omega = 0.0;


	FILE *fp;

	//fp = fopen("Details.dat","at");

	double OldN = N[1][1]; double OldZ = Z[1][1]; double OldD = D[1][1]; 

	int Count = 1;

	/*
	double pt1 = ran2(&seed);
	double pt2 = ran2(&seed);
	double gauss = sigma*(sqrt(-2*log(pt1))*cos(2.0*PI*pt2));
        */
	//double WhiteNoise = gauss; //gauss;

        double WhiteNoise = gsl_ran_gaussian(r2,sigma);
   	for(i=1;i<=MaxRow;i++)
				for(j=1;j<=MaxCol;j++){
					if(N[i][j]>0){
					  
					       if(ForestType[i][j]==0){ //This appears to be non-inducible
						 //Pars[1] *= 1.75;
						 alpha = 1.0; //1.0 is canonical //1.25; //DemogParams[5]*0.01;
						 beta = DemogParams[8];
						  //D[i][j] = 0.0;
						 lambda = DemogParams[0];
					       }else { 
						 //Pars[1] = nubar;
						 alpha = DemogParams[5];
						 lambda = DemogParams[0];
						 beta = DemogParams[8];
					       }
					  	


					        predation = 1 - (2*a*b*N[i][j]/(b*b + N[i][j]*N[i][j]));
						if(Z[i][j]>0){
							StateVars[0] = N[i][j]; StateVars[1] = Z[i][j];
							N2 = N[i][j]; Z2 = Z[i][j]; D2 = D[i][j];
						
								k = EpiParams[2]*exp(theta*(D2+D0));
							
							
							//FractInf = zbrentd2(fx,k,N2,Z2,cx,dx,tol);
							Covars[1] = N2; 
							Covars[2]  = Z2; 
							Pars[2] = k;
							Pars[1] = nubar*exp(-omega*(D2+D0)); //This is to allow defenses to affect nubar, but omega = 0 above is the default

			        			FractInf = Distd1(Pars,Covars);
							if(FractInf<0){ Stop = 1;}
							
						        
							Survivors = (1.0-FractInf)*N2;
	                                           //predation = 1 - (2*a*b*Survivors/(b*b + Survivors*Survivors));
							//Z[i][j] = phi*(N[i][j] - Survivors) + gamma*Z[i][j];
							Z[i][j] = phi*FractInf*N2 + gamma*Z2;
							N[i][j] = lambda*Survivors*predation;
                                                 if(sigma>0)
							  N[i][j] *= exp(WhiteNoise);
							D[i][j] = D2*N2*alpha/(1 + beta*D2);
							//D[i][j] = 0;

							//junk = Z[i][j]; junk = N[i][j];
						} else {
							Z[i][j] = 0;
							N[i][j] = lambda*N2*predation;
                            			if(sigma>0)
							  N[i][j] *= exp(WhiteNoise);
							D[i][j] = D2*N2*alpha/(1 + beta*D2);
						}
					} else {
						Z[i][j] = N[i][j] = 0; D[i][j] = 0;
					}
					Count++;
					
					
					Pop = ((i-1)*MaxRow) + j;
					//N[i][j] *= exp(GaussianNoise[Pop]);

					
				}

				predation = 0;
				

				return;



				
}  //EpiDemog

#include "time.h"

//%SUMATION FUNCTION%
double SumPop(double **StateVar, int MaxRow, int MaxCol){

	int i,j;
	double Sum = 0;
	for(i=1;i<=MaxRow;i++)
		for(j=1;j<=MaxCol;j++)
			Sum += StateVar[i][j];

	return Sum;

}

//%STOCHASTICITY%
int GetNoise(double *GaussianNoise, double sigma, int MaxRow, int MaxCol, long seed){  //this function is pointless, because the noise is site-specific, and is calculated in EpiDemog

	int i,j;
	int Pop;
	long idum;

	
	srand((unsigned)time(NULL)); // initialize ran # generator
 	idum= -rand();
	seed = idum;

	for(i=1;i<=MaxRow;i++)
		for(j=1;j<=MaxCol;j++){
			Pop = ((i-1)*MaxRow) + j;
			//GaussianNoise[Pop] = sigma*gasdev(&sigma);
			GaussianNoise[Pop] = 0.0; //gsl_ran_gaussian(r2,sigma);
		}

	return 0;

}

double Corrln(double ***NStor, int MaxRow, int MaxCol, int tStorMax, float RowScale, float ColScale, double *BinDist, double *CorrlnStor){

	int i,j,k,l;
	int Bin;
	double r;
	double Num;
	long NumCorrln[500];
	FILE *fp;
	double NumBins = 100;

    double Dist, MaxDist, MinDist, jump;

	MaxDist = sqrt(  pow(((1-MaxRow)*RowScale),2) + pow(((1-MaxCol)*ColScale),2)  );

	if(RowScale>ColScale)
		MinDist = RowScale;
	else
		MinDist = ColScale;

	jump = (log(MaxDist) - log(MinDist))/NumBins;
	Bin = 1;
	Dist = log(MinDist);

	while((Dist<log(MaxDist))&&(Bin<=NumBins)){

		BinDist[Bin] = Dist;
		NumCorrln[Bin] = 0;
		CorrlnStor[Bin] = 0;

		Bin++;
		Dist += jump;
		
	}

	Num = 1;
	//fp = fopen("rVals.dat","wt");
	for(i=1;i<=MaxRow;i++)
		for(j=1;j<=MaxCol;j++)
			for(k=i;k<=MaxRow;k++){
				if(k==i) 
					l = k + 1;
				else
					l = 1;
				while(l<=MaxCol){
					pearsn(NStor[i][j],NStor[k][l],tStorMax,&r);

					Dist = sqrt(  pow(((i-k)*RowScale),2) + pow(((j-l)*ColScale),2)  );
					Bin = 1;
					while((log(Dist)>BinDist[Bin])&&(Bin<99))
						Bin++;

					CorrlnStor[Bin] += r;
					//fprintf(fp,"%d\t%f\n",Bin,r);
					NumCorrln[Bin]++;

					

					l++;
				} //while l
			} //k

			//fclose(fp);
			//exit(1);



	//fp = fopen("CorrlnOut.dat","wt");
	for(Bin=1;Bin<=100;Bin++){
		if(NumCorrln[Bin]>0){
			CorrlnStor[Bin] /= NumCorrln[Bin];
			//fprintf(fp,"%d\t%d\t%f\t%f\n",Bin,NumCorrln[Bin],exp(BinDist[Bin]),CorrlnStor[Bin]);
			//fprintf(fp,"%f\t%f\n",exp(BinDist[Bin]),CorrlnStor[Bin]);
			//fprintf(fp,"%d\t%f\n",Bin,CorrlnStor[Bin]);
		}
		
	}
	//fclose(fp);

	return 0;

} //Corrln

//%CALCULATE INVASION SPEED%
double WaveSpeed(double **N, float THold, float OldPoint, float *NewPoint, int MaxRow, int MaxCol){

	// GREG'S ORIGINAL CODE:
	//int j;
	//int StartPoint = ((float) MaxRow)/2.0; //Greg's; don't know why this doesn't match values in Initialize
	//double junk;

	//j = MaxCol;
	//do{ 
	//	j--; 
	//	junk = N[j][StartPoint];	    
	//}while((N[j][StartPoint]<2.0)&&(j>1));

	//*NewPoint = j;
	
	//return *NewPoint - OldPoint;		
	//END GREG'S CODE
	
	int startRow = 1;
	int startCol = 1;
	int i, j;
	double howfar;
	double junk = 0;

	for(i = 1; i < MaxRow+1; i++){
		for(j = 1; j < MaxCol+1; j++){
			if (N[i][j] > THold){
				howfar = sqrt(  pow((double)(i-startRow),2) + pow((double)(j-startCol),2)  );
				if (howfar > junk){
					junk = howfar;
				}
			}
		}
	}

	*NewPoint = junk;

	return *NewPoint;
	

}
		

#include "time.h"
long idum;


int main(void)
{

	
	int i,j;
	FILE *fp, *fp2;
	int GridSize = 100; //100;//100 goes with J //50; //grid has GridSize X GridSize cell
	int MaxRow = GridSize; int MaxCol = GridSize; //dimensions of the spatial grid
	float RowDim = 25; //30; //1000; 50 if grid is 25, otherwise 100 if grid is 100  //horizontal scale (total number of km in GridSize (or btwn 1 and MaxRow+1))
	float ColDim = RowDim; //vertical scale
	float RowScale = RowDim/((float) MaxRow);  //kilometers per grid point, horizontal
	float ColScale = ColDim/((float) MaxCol); //kilometers per grid point, vertical

	//demographic parameters:
	double kernelCOEF = 18.916168; //dispersal parameter for normal kernel (m^-1) for ballooning larvae
	double LongDistParam = 1/0.175;//from Pauline's dissertation
	double LongDistFract = 0.00001;//from Pauline's dissertation
	double truncateProb = 1e-10;  //max dispersal allowed is the distance with probability of truncateProb
	double t, maxt;
	int tStor, tStorMax = 250;
	int Rlzn, MaxRlzns;
	float kParam = 1.06;
	double InitN = 1; 
	double InitZ = 1;
	double InitD = 1;
	float InitFreq = 0.5;
	float phi = 3, lambda = 70.0;
	float a = 0.97;
	float b = 0.14;
	float gamma = 0.0;
	double sigma; //1.5;//tuned to match synchrony data in Peltonen et al.
	float CorrlnDistParam = 1e-5; //3.3e-3;//tuned to match synchrony data in Peltonen et al.	
	double rho = 1;
	float alpha, theta, D0, beta;

	
//lam <- 5.0; phi <- 1.125; gam <- 0.125; V <- 200; k <- 1/V 
//alpha <- 20; Beta = 15; Co = 5.0; omega <- 0.0;  psi <- 0.00; xi = 0.00;
	lambda = 74.6; phi = 0.5; gamma = 0.3; 
	kParam = 1.0/60.0; 
	alpha = 20; 
	D0 = 5; 
	//InitD = D0;
	beta = 15; 
	theta = 1;
	a = 0.967; b = 0.14; 
 
//from SOM 21 March 2010
	
	theta = 1;
	lambda = 5.1; phi = 1.125; gamma = 0.225;
	alpha = 20; beta = 15; 
	a = 0.775; b = 0.14;
	//V = 200; 
	D0 = 4.95; //5
	//InitD = D0;
	kParam = 1.0/65.0; //1.0/65.0;
	double d = 0.00;

	
	//From Bret, 12 May 2010
lambda  = 5.5; phi = 1.125; gamma = 0.225; kParam = 1.0/52.0;
	

InitD = 0.01; InitN = 0.1; InitZ = 0.1; 
 //printf("InitN:%e InitZ:%e InitD:%e\n",InitN,InitZ,InitD);



//Testing - current parameters, really
//From R code:
//Old Set 1: alpha 4, beta 100, phi 2.5, kP 0.04, gamma 0.1, D0 3.1, sigma 0.025;
//Old Set 2: alpha 5, beta 100, phi 1.5, kP 0.04, gamma 0.2, D0 3.2, sigma 0.025
 

 lambda = 74.6; a = 0.95; b = 0.14; 
 alpha = 2.5; //8.0; //8 is the new canonical //4 is canonical //5; //4 is default, 3.5 works for X200
 beta = 100.0; 
 phi = 0.5; //0.5 is what  you've been using //1.5; //1; //2.5; 
 kParam = 0.04; //0.04; is canonical
 gamma = 0.2;  
 D0 = 3.0; //3.2; //3.1; //InitD was 4.25
 //THE TRUE RE-SCALING IS AS FOLLOWS: (this comment appears to be nonsense)

 sigma = 0.3; //0.2; //0.3; is canonical

 theta = 1; 
 MaxRlzns = 1;
 maxt = 100;
 float IndProb = 0.3;
 int Step;
 if(MaxRlzns>1){
	 Step = 1;
 }else{
	 Step = 10;
 }
 int TS = 0;

 //Gives a superharmonic: lambda 74.6, a 0.9, b 0.14, alpha 4, beta 80, gamma 0.2, D0 4.0, sigma 0, theta 1, phi 5
 //Lower period: phi 3, 0.25
 //lower period, multiple harmonics: phi 2, alpha 4, D0 4.1, nice super phi 2, alpha 2, D0 4.1

 //Gives period 9: lambda 50, alpha 2, beta 16.2, phi 8, InitD 4.35, sigma 0.01
 
//A - alpha 200, beta 200.  B is A, only Prob = 0.5; C has alpha 2e-4, beta 2e-4, Prob = 0; D has alpha 200, beta 200, Prob = 0.75; E is D with Prob = 0.6
//F - is with predation AFTER infection, alpha 0.305, beta 0.012; G alpha 2, beta 0.1; H alpha 10, beta 10; I alpha 2e-3, beta 2e-3; J alpha 2e-2 beta 2e-3; L alpha 0.5, beta 5e-3
//M alpha 2e-3, beta 2e3; N has alpha 2e-4, beta 2e3, P is like N only phi 1 not 1.5; Q like P, only alpha 2e2, beta 2e2; R phi 1, alpha 2e3, beta 2e2

//H series is doing a grid of values of alpha and beta
//1 has phi 1, Prob 0. 
//2A has phi 10, prob 0
//3A has phi 50, prob 0; B has phi 50, prob 0.75; C has phi 50, prob 0.95; D is sort of
//E is on a smaller scale (25 km) and lower stoch (sigma = 0.1); F has Long dist disp'l turned off; G like E only sigma is 1.5, to enhance synch
//H is like D, only sigma = 0.25; J like H, only Rows=Cols = 100
//K is 25 X 25

//TS is 0.6, TSB is 0.0, TSC is 0.9, TSD is 0.2, TSE is 0.8, TSF is 0.1, TSG is 0.01, TSH is 0.05, TSJ is 0.01 with sigma = 0.3
//Low is 0.25, Low2 is 0.5, Low3 is 0.75, times either alpha or lambda
//Specific is alf 1 lowlam 0.75, B is alf 0.5

//D is alf 1, lowlam 1, E is alf 1, lowlam 0.75, 50 rlzns
//Specific means only particular percentages of oak, rather than a whole range.  46 means only 46%, 21 means only 21%
//Sp is short for "Specific", G means following numbers are grid sizes, these are for 1600 (40X40), and for 2500 (50X50)
//std is with RowDim 30, B is with RowDim 60, C like std, but lowlam = 1.0 (and lowalf 1), D again has default lowlam (0.75), but lowalf = 0.5; 
//E is D with phi = 0.25 instead of 0.5 (D, E 25 rlzns, others50)
//F has Rowdim 30, lowlam = 0.75, lowalf = 1, phi = 0.25, G like F, but lowlam = 1.0, lowalf = 0.5, phi = 0.5; H like G, but lowalf = 1.25; I like H, but lowlam = 0.75
//J like C, only lowalf 1.25, K like C, only lowalf 1.5; L like K, only lowlam 0.75; M like L, only lowalf 2, lowlam 0.75; N like M, only lowlam 1
//P is C but with Rowdim 40, Q has Rodim 20; R has RowDim 25; S like R, but lowalf 0.5; T is 21 R, which got screwed up (IndProb 0.81 instead of 0.79)
//21 is prob of oak, so IndProb is 0.81 (NO! 0.79), 46 corresponds to IndProb 0.54
//U has RowDim 25, lowalf 1, lowlam 0.75; V has alf 4.5 (NOT lowalf) lowlam 1, RowDim 25; W has alf 4, lowlam 1, RowDim 25, sig 0.3
//X has RowDim 25, lowalf 1, lowlam , but alf 3.5 (U and V 6e4 are almost certainly 2e4); Y like X, only alpha = 3.75
//X125 means 125 X 125, as opposed to G1e4, which means 100 X 100; Z like W, only sig 0.4; A like Z, only sig 0.2
 //W? and Z had alpha only 3.75
 //A has sig 0.3, alpha 4.0; B has sig 0.3, alpha 3.5; C has sig 0.3, alpha 3.25; D like B but sig 0.4
//LoLA is low lambda A, which low lambda.  A is 0.75X in non-inducible, B is 0.5, C is 0.25, all with R params,
//NoAlf means no alpha diffs, but low lambda.  A is 075X, B is 0.5, C is 0.25, all with R params; 
//For 15/43, A is actually 0.25 (whoops), and C is 0.75; D is 0.05; 
//X has low alpha at 0.5, lowlam 0.75, XB has low alpha 1.5, lowlam 0.75; XC has low alpha 1, lowlam 0.75, hi alpha 5; XD is C but hi alpha 5.5: E has hi alpha 4, lowlam 0.85
//F has has hi alpha 3.5, lowlam 0.85; G has hi alpha 4, lowlam 0.9; H has low beta 0.2X, lowlam 0.9; I has low beta 0.5X, lowlam 0.9, J has lowbeta 2.5X, lowlam 0.9
//XK has low beta 5.0X, lowlam 0.9; L has low beta 1X, lowlam 0.75, low alpha 0.2; M has low beta 1X, lowlam 0.75, low alpha 0; XN is low alpha 0, lowlam 0.75, high alpha 2.0
//XC43 is actually XC15!!
//NoD is no dispersal
//TURN DISPERSAL BACK ON!  FOR GOD'S SAKE!!!!!!!!!!!!!!!!

//Fix:
//Original is SpR43G1e4, except fract oaks is 0.34.  A is SpR43G1e4 with IndProb really 0.43.  B has alpha 5; C has hi alpha 4, low alpha 0.5
//D has alpha 4, low alpha 1, phi 1.0 instead of 0.5; E is D with gamma 25/12 instead of 15/12; F is D with gamma 56/10 instead of 15/12; G like F with phi 0.75 instead of 0.5
//H is E with alpha = 5; I has k = 0.03, alpha = 5; J has two nubars, non-inducible is 2X; K like J, only no change in k
//L has low nubar 0.5X, high alpha 2.5; M is L with high alpha 1.5; N is is L with high alpha 1, low alpha 0.25; 
//O is E with high alpha 4, low alpha 1, low nubar 1.05*high; P is O with low nubar 1.25*high; Q is P with low alpha 1.5; R is P with low alpha 2.0; 
//S is O with low nubar 1.3*high, low alpha 1.5; T is S with low nubar 1.4*high; U is T with low alpha 1.25; V is U with phi = 0.25 to reduce crashing...
//W has high nubar 1.75X, high alpha 1.1, low alpha 1.0
//EB is a repeat of E, but phi 0.5; EC is EB with alpha 4.5; ED is EB with alpha 4, phi 1 (sigma 0.3 in all cases); EF (or EFG3 instead of EFG43) with alpha 5, phi 0.5; EG alpha 5.5, phi 0.5
//EH has alpha 6, phi 0.5, sigma 0.3; EI has alpha 6.5; EJ alpha 7; EK alpha 8; EL alpha 9; EM alpha 11

//NoLg: prints out N, not log10(N)
//EC is FixEC (alpha 4.5); EN alpha 13; EI is FixEI, alpha 6.5; EK akpha 8; EO is EK with time step h = 0.002; EP has hi alpha 8, lo alpha 0.8; EQ like EP hi alpha 10; ER like EP hi alpha 12
//ES has low alpha 0.5, hi alpha 10; ET has low alpha 0.25, hi alpha 10
//X200: EL has hi alpha 6, EM hi alpha 10 both have low alpha 1, EN has alpha 4 
//EKB is a repeat of EK

//nu2 has nubar = exp(-omega*(D2+D0));
//B has omega 0.2 alpha 12 sigma 0.3, C has omega 0.2 alpha 10 sigma 0.2; D has omega 0.17 alpha 10 sigma 0.2; E has omega 0.17 alpha 8 sigma 0.2; F has omega 0.17 alpha 9 sigma 0.2; G has omega 0.17 alpha 12 sigma 0.2
//H like G with alpha 6; J like G with alpha 14; K like G with alpha 16; L has omega 0.05, alpha 8; M has omega 0.1, alpha 8; N has omega 0.1, alpha 10; O has omega 0.1, alpha 12; P has omega 0.07, alpha 10
//Q  has omega 0.1, alpha 9; R has omega 0.1, alpha 8.5 (always sigma = 0.2)

//LoL2: lower lambda in non-induced plots, alpha 8, omega 0, sigma 0.3
//A has low alpha 1.0, hi alpha 8, low lambda 0.95X; B like A with low lambda 0.9; C like A with low lambda 0.85; D like A with low lambda 0.8; E like A with low lambda 0.75; F like A with low lambda 0.7
//G has low lambda 0.75, hi alpha 11; H has low lambda 0.75, hi alpha 13; I has low lambda 0.7, hi alpha 13; J has low lambda 0.7, hi alpha 14; K has low lambda 0.65, hi alpha 13; L like K with hi alpha 14
//M like K with hi alpha 14, low alpha 1.5; N like K with hi alpha 16, low alpha 1.5; O with low lambda 0.65, hi alpha 16, low alpha 1; P like O, hi alpha 18
//EA has hi alpha 9, low alpha 1, low lambda 0.75 (like E, with alpha 9); EB like EA with low lambda 0.7, high alpha 11; GB like G with alpha 12; GC like G with alpha 13; GD alpha 10; GE alpha 9
//GF like G with alpha 8; GG alpha 7; GH alpha 6
//DA has low lambda 0.8X, high alpha 9; DB has alpha 10; DC has alpha 11; DD has alpha 12; DE has alpha 7; DF has alpha 6; DG has alpha 5
//Q has low alpha 0.5, high alpha 8, low lambda 0.85;

//nuLoL: BOTH variable nu and low lambda on non-inducible plots
//A has omega 0.8, alpha 10, low lambda 0.75X

//NoAlf2: alpha low same as alpha high, but with lower lambda
//A has alpha 8, lo lambda 0.9; B has alpha 1, lo lambda 0.9; C like B, lo lambda 0.8; D like C, lo lambda 0.7

//NLG - same as NoLg, except uses pid's instead of time/date stamp
//EKC and EKD like EK
//LoAlf25 - alpha on non-induced is 0.25 instead of 1, 50 is alpha 0.5
//HiAlf4 - alpha on induced is 4 instead of 8, HiAlf25 is actually alpha 2.5



	char *test[] =
	  {
	    "HiAlf25EKD43X100", //"NLGEKD43X100", //"NoAlf2C43G1e4", //"LoL2Q43G1e4", //"nu2R43G1e4", //"NoLgEK43G1e4", //"SpRC43G1e4", //0
	    "NDNHiAlf25EKD43X100", //"NDNLEKD43X100", //"NDNoAlf2C43G1e4", //"NDLoL2Q43G1e4", //"NDnu2R43G1e4", //"NDNoLgEK43G1e4", //"NDTotSpRC43G1e4", //1
	    "DNHiAlf25EKD43X100", //"DNLEKD43X100", //"DNoAlf2C43G1e4", //"DLoL2Q43G1e4", //"Dnu2R43G1e4", //"DNoLgEK43G1e4", //"DRCSp43G1e4", //2 Defense
	    "FTHiAlf25EKD43X100" //"FTNLEKD43X100" //3
	  };

char testbuff2[128];
char buffer0[128],buffer1[128],buffer2[128],buffer3[128];
char bufferB[128];
char fname;
int pid;
pid = getpid();
//printf("pid:%d\n",pid);




	char *strFileType = ".dat";
	strcpy(buffer0, test[0]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer0, bufferB);
	strcat(buffer0,strFileType); 
	test[0] =  buffer0;

	strcpy(buffer1, test[1]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer1, bufferB);
	strcat(buffer1,strFileType); 
	test[1] =  buffer1;

	strcpy(buffer2, test[2]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer2, bufferB);
	strcat(buffer2,strFileType); 
	test[2] =  buffer2;

	strcpy(buffer3, test[3]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer3, bufferB);
	strcat(buffer3,strFileType); 
	test[3] =  buffer3;

	
/*
	//for(i=0;i<=3;i++){
     		strcpy(buffer0, test[0]);
     		sprintf(buffer2, "%d", pid);
     		strcat(buffer0, buffer2);
		strcat(buffer0,strFileType); 
		test[0] =  buffer;

		
		//strcpy(testbuff2,test[i]);
		
		
	//}
*/


/*
	for(i=0;i<=3;i++){
		strcpy(testbuff2,test[i]);
		fp = fopen(testbuff2,"a");
		fprintf(fp,"Testing2...\n");
		fclose(fp);	
	}
exit(1);
*/
		


//exit(1);
			

		
/* This is old code that slaps time/date on end of output files, now superseded (I hope) by code that slaps pid on end instead
	//FILE *fp1;
	struct tm *ptr; //struct defined in time.h;  used to hold the time and date
	time_t lt; //type defined in time.h; for storing the calendar time.
	char strDate[30]; //string to hold date and time
	char *strFileName = "file_name_"; 
	char *strFileNameDate; 
	
	//int i;

	//wait(60);
		lt = time(NULL);
		ptr = localtime(&lt); //converts calendar time to local time

		strftime(strDate, 30, "%M%S", ptr); //adds to strDate month, day, year, hour, minute, and second from ptr
		//strftime(strDate, 30, "%m_%d_%y_%H:%M:%S", ptr); //adds to strDate month, day, year, hour, minute, and second from ptr
		//strftime(strDate, 30, "%S", ptr); //adds to strDate month, day, year, hour, minute, and second from ptr
                //strftime(strDate, 30, "%S_%M_%H",ptr);

	// printf(str1);
	for(i=0;i<=3;i++){
	 
		strFileNameDate = (char*)calloc((strlen(test[i])+ strlen(strFileType)+strlen(strDate)+1), sizeof(char)); 
	//allocate memory big enough to hold all 3 strings

		strcat(strFileNameDate,test[i]); //add strFileName to 1st empty memory space in strFileNameDate
		strcat(strFileNameDate,strDate); //add strDate to 1st empty memory space in strFileNameDate
		strcat(strFileNameDate,strFileType); 
	
		//printf("\n file name = %s\n date = %s \n file type = %s \n file name with date and type =  %s\n", strFileName, strDate, strFileType, strFileNameDate);

		test[i] = strFileNameDate;
	}
*/

	
	/*
	fp2 = fopen("TwoDimSpace.dat","wt");
	fclose(fp2);
	*/

	fp = fopen("TotPop.dat","w");
	fclose(fp);

	//fp = fopen("SpectraOut20A.dat","w");
	//fclose(fp);

	//fp = fopen("junk.dat","w");
	//fclose(fp);

	/*
	fp = fopen("Noise.dat","wt");
	fclose(fp);

	fp = fopen("NoiseII.dat","wt");
	fclose(fp);
	*/

	char testbuff[128];
	//A has alpha 200, beta 200; note that for G20 (20 X 20), there is no A	
        //B has phi 1.5, C has alpha 20, beta 20 (C apparently has B in its file names, D is like B, only Dim = 600, not 100; E is D with fewer rlzns
        //30A has phi 1.5, alpha 200, beta 200, Dim 100

		strcpy(testbuff,test[0]);
			fp = fopen(testbuff,"w");
			fclose(fp);


			/*
		strcpy(testbuff,test[1]);
			fp = fopen(testbuff,"w");
			fclose(fp);
			*/

	strcpy(testbuff,test[2]);
			fp = fopen(testbuff,"w");
			fclose(fp);

			/* Not in use
	strcpy(testbuff,test[3]);
			fp = fopen(testbuff,"w");
			fclose(fp);
			*/


	//demographic parameters end	
	
	//long seed = -2;

	float junk;

	//N is host pop size, Z is pathogen pop size, D is the defense
	double **N;
	double **Z;
	double **D;
	double **ForestType;
	double ***NStor;
	double **NStorOneD;
	
	//variables for storing parameters and distances and so forth
	//double Dist;
	//double MovingN, MovingZ;
	float *EpiParams, *DemogParams;
	double *StateVars;
	//double Survivors;


	//storage parameters for storing data
	double SumN[1000];
	double SumZ[1000];
	double SumD[1000];
	double *GaussianNoise;
	//int Pop;
	int NumPops = MaxRow*MaxCol;
	double *BinDist, *CorrlnStor;
	double *AvgCorrln, *StdDevCorrln;	

	double **Sigma;
	int number = MaxRow*MaxCol;

	//Masting parameters //(not used)
	//float Mean = 4;
	//int NextMastYr = poidev(Mean,&seed);
	//int MastYrCount = 0;
	//float atemp;

	//Speed parameters
	float THold = b;//b is the "Allee" threshold for the predation model
	float OldPoint, NewPoint;
	float Speed;
	int StartPoint;
	double *GaussRandNums;
	double pt1,pt2;

	//seed = SeedGet();

	srand((unsigned)time(NULL)); // initialize ran # generator
 	idum= -rand();
	long seed;
	seed = idum;
       seed = pid; //Why not?
	
	//seed = -2;

	//Allocating Memory below//
	
	GaussRandNums = dvector(1,NumPops);
	StateVars = dvector(0,10);
	EpiParams = vector(0,10); DemogParams = vector(0,20);
	GaussianNoise = dvector(0,NumPops);

	N = dmatrix(0,MaxRow+1,0,MaxCol+1);
	Z = dmatrix(0,MaxRow+1,0,MaxCol+1);
	D = dmatrix(0,MaxRow+1,0,MaxCol+1);
	
	ForestType = dmatrix(0,MaxRow+1,0,MaxCol+1);
	NStor = d3tensor(0,MaxRow+1,0,MaxCol+1,0,tStorMax+1);
	NStorOneD = dmatrix(0,maxt,0,MaxCol+1);

	CorrlnStor = dvector(0,500);
	AvgCorrln = dvector(0,500);
	StdDevCorrln = dvector(0,500);
	BinDist = dvector(0,500);

	Sigma = dmatrix(0,number,0,number);

	//Parameter storage
	EpiParams[0] = 0; EpiParams[1] = 0; 
	EpiParams[2] = kParam;
	DemogParams[0] = lambda; DemogParams[1] = phi; DemogParams[2] = gamma;
	DemogParams[3] = a; DemogParams[4] = b;
	DemogParams[5] = alpha; DemogParams[6] = theta; DemogParams[7] = D0; DemogParams[8] = beta;
        DemogParams[9] = sigma;

	double sigma2 = 0.0; 
	//CorrlnMatrix(rho,0.0,sigma2,seed,MaxRow,MaxCol,RowScale,ColScale,CorrlnDistParam,Sigma); //Don't need this, right?

	//fp = fopen("WaveSpeed.dat","wt");
	//fclose(fp);

	for(i=1;i<=100;i++){
		AvgCorrln[i] = 0.0;
		StdDevCorrln[i] = 0.0;
	}

			
 const gsl_rng_type *T2;
 gsl_rng_env_setup();
 srand((unsigned) time(NULL));
 seed = -rand();
 seed = pid;
 gsl_rng_default_seed = seed;

 T2 = gsl_rng_default;
 r2 = gsl_rng_alloc(T2);

 fp = fopen("maxtB.dat","w");
 fprintf(fp,"%f\n",maxt);
 fclose(fp);
 
 //fp = fopen("H3K.dat","w");
 //fclose(fp);

 //fp = fopen("Details.dat","w");
 //fclose(fp);

 //for(alpha=20.0;alpha<=35;alpha+=1.0e10)
 //for(beta=35;beta<=1e2;beta*=1.0e10)
 for(IndProb=0.57;IndProb<=0.795e1;IndProb+=0.25e10){ //REally it's 0.54 or 0.81 (NO! it's 0.57 or 0.85 as of March 2012), dunno about this other stuff->Or 0.54, 0.795, +=0.25; I think that the idea is that IndProb is actually prob of NON induction
	for(Rlzn=1;Rlzn<=MaxRlzns;Rlzn++){

	  //fp = fopen("AlphaBeta.dat","a");
	  //fprintf(fp,"%f  %f\n",alpha,beta);
	  //fclose(fp);



		double Lower = 0.1; double Upper = 2.0;
		InitN = Lower + (Upper-Lower)*gsl_rng_uniform(r2);
		InitZ = Lower + (Upper-Lower)*gsl_rng_uniform(r2);
		//InitD = Lower + (D0-Lower)*gsl_rng_uniform(r2);
		

		Initialize(N,Z,D,ForestType,InitN,InitZ,InitD,MaxRow,MaxCol,IndProb);  //initializing pop densities
		double SumBalloonFract = InitBalloon(kernelCOEF, MaxRow, MaxCol, RowScale, ColScale);
		double SumLDDFract = InitLDD(LongDistParam, LongDistFract, MaxRow, MaxCol, RowScale, ColScale, truncateProb);

		//Writing ForestType to a file
		strcpy(testbuff,test[3]);
		fp = fopen(testbuff,"a");
		fprintf(fp,"%f  %f  ",alpha,IndProb);
		for(i=1;i<=MaxRow;i++)
			for(j=1;j<=MaxCol;j++){
			    fprintf(fp,"%f ",(ForestType[i][j]));
			}
		fprintf(fp,"\n");
		fclose(fp);
		
		tStor = 1;
		for(t=1;t<=maxt;t++){

			//printf("Rlzn:%d t:%f\n",Rlzn,t);

	  /*
			fp = fopen("TotPop.dat","a");
			fprintf(fp,"%f %e %e %e\n",t,log10(N[1][1]),log10(Z[1][1]),log10(D[1][1]));
			fclose(fp);
		  */

		  //if(fmod(t,25)==0){
		  //fp = fopen("RlznCountB.dat","a");
		  //fprintf(fp,"%d %f\n",Rlzn,t);
		  //fclose(fp);
		    //}
		  /* This is really for long time-series only
		  	strcpy(testbuff,test[1]);
			fp = fopen(testbuff,"a");
			fprintf(fp,"%f %f %f %f\n",t,(SumPop(N,MaxRow,MaxCol)),(SumPop(Z,MaxRow,MaxCol)),(SumPop(D,MaxRow,MaxCol)));
			fclose(fp);
		  */

		  /*
		  fp = fopen("Noise.dat","a");
		  fprintf(fp,"Rlzn:%d t:%f \n",Rlzn,t);
		  fclose(fp);
		  */
	
		   //Next 2 lines are dispersal.  Un-comment them to turn dispersal back on
		  DispersalFunc2(N, Z, SumLDDFract, LongDistParam, LongDistFract, MaxRow, MaxCol, RowScale, ColScale, truncateProb);  //LDD
		  DispersalFunc(N, Z, SumBalloonFract, kernelCOEF, LongDistFract, MaxRow, MaxCol, RowScale, ColScale);  //ballooning
	
		//DispersalFuncOttar(N, Z, SumBalloonFract, kernelCOEF, d, MaxRow, MaxCol, RowScale, ColScale);  //Don't do this one!  It's Ottar's dispersal function
			
			for(j=1;j<=NumPops;j++){

			  //pt1 = ran2(&seed);

			  //pt2 = ran2(&seed);
			  double gauss = 	0.0; // gsl_ran_gaussian(r2,sigma); //Stochasticity is actually ELSEWHERE, specifically inside EpiDemog.  This chunk appears to be pointless.
				GaussRandNums[j] = gauss;
				  //(sqrt(-2*log(pt1))*cos(2.0*PI*pt2));
				//GaussRandNums[j] = 0.0;
				if(sigma<=0) GaussianNoise[j] = 0.0; 

			}
			
			
			/*
			if(sigma>0){
				DistCorrldRand2(rho,0.0,sigma,seed,MaxRow,MaxCol,RowScale,ColScale,CorrlnDistParam,Sigma,GaussianNoise,GaussRandNums);
			}
			*/
			
		
			
			DemogParams[5] = alpha; DemogParams[8] = beta;
                        //printf("gam:%f lambda:%f a:%f b:%f kParam:%f D0:%f alpha:%f beta:%f \n",gamma,lambda,a,b,kParam,D0,alpha,beta);
	                //printf("Before   N:%e Z:%e D:%e\n",N[1][1],Z[1][1],D[1][1]);
			EpiDemog(N,Z,D,ForestType,MaxRow,MaxCol,DemogParams,EpiParams,GaussianNoise);  //the epidemic




			//if(fmod((double) t,Step) == 0) 
			  //printf("IndProb:%f Rlzns:%d Gen:%f SumN:%e SumZ:%e SumD:%e\n",IndProb,Rlzn,t,(SumPop(N,MaxRow,MaxCol)),(SumPop(Z,MaxRow,MaxCol)),(SumPop(D,MaxRow,MaxCol)));

			

			
			//fp = fopen("OneDimSpace.dat","at");
			//for(j=1;j<=MaxCol;j++){
			//StartPoint = ((float) MaxCol)/2.0;
			//if(N[j][StartPoint]>0)
			//NStorOneD[((int) t)][j] = log10(N[j][StartPoint] + 1);
			//else 
			//NStorOneD[((int) t)][j] = 0;
			//}
			//fclose(fp);


			if((maxt-t)<=tStorMax){

				for(i=1;i<=MaxRow;i++)
					for(j=1;j<=MaxCol;j++)
						NStor[i][j][tStor] = log10(N[i][j]);
				tStor++;
			}

			//Next loop is for storing 2-D space
			/*
			if(t>0){
			
			//strcpy(testbuff,test[1]);
			//fp = fopen(testbuff,"a");
			
			fp = fopen("TwoDimSpace.dat","at");
			//fp = fopen("junk.dat","at");
			
			for(i=1;i<=MaxRow;i++){
				for(j=1;j<=MaxCol;j++){
					//fprintf(fp2,"%d\t",N[i][j]);
					if(N[i][j]==0)
						//fprintf(fp2,"NA\t");
						fprintf(fp,"0\t");
					else{
						//fprintf(fp,"%f\t", log10((double) N[i][j])  );
						fprintf(fp,"%f\t", ((double) N[i][j])  );
					}
					//fprintf(fp,"%f %f %f\n",N[i][j],Z[i][j],D[i][j]);
					//getc(stdin);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
			
			//fclose(fp);
			} //t > 0
			


			Speed = WaveSpeed(N,THold,OldPoint,&NewPoint,MaxRow,MaxCol);
			*/
			
			//OldPoint = NewPoint;

			//fp2 = fopen("WaveSpeed.dat","at");
			//fprintf(fp2,"%f\t%f\t%f\n",t,Speed*ColScale,(NewPoint-1)*ColScale);
			//fprintf(fp2,"%f\t%f\n",t,Speed*RowScale);
			//fclose(fp2);

		       if(t>TS){
			strcpy(testbuff,test[0]);
			fp = fopen(testbuff,"a");
			//printf("alpha:%f beta:%f\n",alpha,beta);
                        fprintf(fp,"%f  %f  ",alpha,IndProb);
			//fp2 = fopen("SpectraOut20A.dat","at");
			double NTot = 0.0;
			for(i=1;i<=MaxRow;i++)
			  for(j=1;j<=MaxCol;j++){
			    fprintf(fp,"%f ",(N[i][j]));
					NTot += N[i][j];
			  }
					
			fprintf(fp,"\n");
			fclose(fp);
                       }

			strcpy(testbuff,test[2]);
			double DTot = 0.0;
			fp = fopen(testbuff,"a");
			//printf("alpha:%f beta:%f\n",alpha,beta);
                        fprintf(fp,"%f  %f  ",alpha,IndProb);
			//fp2 = fopen("SpectraOut20A.dat","at");
			for(i=1;i<=MaxRow;i++)
			  for(j=1;j<=MaxCol;j++){
			    fprintf(fp,"%f ",D[i][j]);
			    DTot += D[i][j];
			  }
			fprintf(fp,"\n");
			fclose(fp);

			/*
			strcpy(testbuff,test[1]);
			fp = fopen(testbuff,"a");
			fprintf(fp,"%f %f %f\n",t,DTot,DTot);
			fclose(fp);
			*/
			
		} //t

		/*
		junk =  Corrln(NStor,MaxRow,MaxCol,tStorMax,RowScale,ColScale,BinDist,CorrlnStor);

		fp = fopen("CorrlnOut2b.dat","at");
		for(i=1;i<=100;i++){
			AvgCorrln[i] += CorrlnStor[i];
			StdDevCorrln[i] += CorrlnStor[i]*CorrlnStor[i];
			if(CorrlnStor[i]!=0)
				fprintf(fp,"%d\t%d\t%f\n",Rlzn,i,CorrlnStor[i]);

		}
		fclose(fp);
		*/

		//if(Rlzn==2) exit(1);
		
	} //Rlzn
 } //IndProb

	/*
	fp = fopen("AvgCorrln.dat","wt");
	for(i=1;i<=100;i++){
		if(AvgCorrln[i]>0){
			AvgCorrln[i] /= (float) MaxRlzns;
			StdDevCorrln[i] /= (float) MaxRlzns;
			StdDevCorrln[i] -=  AvgCorrln[i]*AvgCorrln[i];
			fprintf(fp,"%f\t%f\t%e\n",exp(BinDist[i]),AvgCorrln[i],StdDevCorrln[i]);
		}
	}
	fclose(fp);
	*/

	

	/*
	
	fp2 = fopen("2DSpace.dat","wt");
	for(i=1;i<=MaxRow;i++)
		for(j=1;j<=MaxCol;j++)
			fprintf(fp2,"%d\t%d\t%f\t%f\n",i,j,N[i][j],Z[i][j]);
	fclose(fp2);


	
	fp = fopen("OneDimSpace.dat","wt");
	for(j=1;j<=MaxCol;j++){
		fprintf(fp,"%d\t",j);
		for(t=1;t<=maxt;t++)
			fprintf(fp,"%f\t",NStorOneD[((int) t)][j]);
			fprintf(fp,"\n");
	}
	fclose(fp);
	*/

	free_dmatrix(N,0,MaxRow+1,0,MaxCol+1);
	free_dmatrix(Z,0,MaxRow+1,0,MaxCol+1);
	free_dmatrix(D,0,MaxRow+1,0,MaxCol+1);
	free_dmatrix(ForestType,0,MaxRow+1,0,MaxCol+1);
	free_d3tensor(NStor,0,MaxRow+1,0,MaxCol+1,0,tStorMax+1);

	free_dvector(StateVars,0,10);
	free_vector(EpiParams,0,10);
	free_vector(DemogParams,0,20);

	free_dmatrix(NStorOneD,0,maxt,0,MaxCol+1);

	free_dvector(BinDist,0,500);
	free_dvector(CorrlnStor,0,500);

	free_dvector(AvgCorrln,0,500);
	free_dvector(StdDevCorrln,0,500);

	free_dmatrix(Sigma,0,number,0,number);

	free_dvector(GaussRandNums,1,MaxRow*MaxCol);
	

	return 0;

	
}

