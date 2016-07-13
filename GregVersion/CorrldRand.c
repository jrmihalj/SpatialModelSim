#define PI 3.1415926535898

void choldc(double **a, int n, double p[])

{

	void nrerror(char error_text[]);

	int i,j,k;

	float sum;

   

	for (i=1;i<=n;i++) {

		for (j=i;j<=n;j++) {

			for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];

			if (i == j) {

				if (sum <= 0.0)

					nrerror("choldc failed");

				p[i]=sqrt(sum);

			} else a[j][i]=sum/p[i];

		}

	}

}

void CorrlnMatrix(float corrln, float mean, float SD, long seed, int MaxRow, int MaxCol, float RowScale, float ColScale, float DistParam,double **Sigma){

	int i,j,k,l;
	int Pop1, Pop2;
	int number = MaxRow*MaxCol;
	float dummy;

	//int number=100; // number of random numbers to generate

	double pt1,pt2,*Z, **C,*Cdiagonal;
	float covar, var;
	float Dist;
	var = SD*SD;

	
	covar = corrln*var;

	for(i=1;i<=MaxRow;i++){

		for(j=1;j<=MaxCol;j++){

			for(k=i;k<=MaxRow;k++){
				/*
				if(k==i) 
					l = k + 1;
				else
				*/
					l = 1;
				while(l<=MaxCol){

				  Dist = sqrt(  pow(((double) (i-k)*RowScale),2) + pow(((double) (j-l)*ColScale),2)  );

					Pop1 = MaxRow*(i-1) + j;
					Pop2 = MaxRow*(k-1) + l;

					dummy = exp(-DistParam*Dist); 
					dummy = Dist;
					dummy = 0;
					if((i==k)&&(j==l))
						Sigma[Pop1][Pop2] = Sigma[Pop2][Pop1] = var;
					else
						Sigma[Pop1][Pop2] = Sigma[Pop2][Pop1] = var*exp(-DistParam*Dist);  //Covariance between random variables

					l++;
				} //l
			}  //k
		} //j
	} //i



	return;

}




void DistCorrldRand(float corrln, float mean, float SD, long seed, int MaxRow, int MaxCol, float RowScale, float ColScale, float DistParam, double **Sigma, double *X){


// Choose Z1,Z2... from N(0,1), Convert to Correl. Normal. random numbers  X1,X2,... 

// by X = C Z  where C is cholesky decomp of the Covariance matrix of X


	/*
	double corrln = 0.999;
	double mean=0, var=0.25, covar;
	long seed = -1;
	*/

	FILE *fp;

	int i,j,k,l;
	int Pop1, Pop2;
	int number = MaxRow*MaxCol;

	//int number=100; // number of random numbers to generate

	double pt1,pt2,*Z, **C,*Cdiagonal;
	float covar, var;
	float Dist;
	long junk;

	var = SD*SD;


	covar = corrln*var;

	//X=dvector(1,number);  // Vector of 10 correlated normal random numbers

	Z=dvector(1,number);   // Vector of 10 normal random numbers

	//Sigma=dmatrix(1,number,1,number);  //Matrix of variances and covariances 

	C=dmatrix(1,number,1,number);

	Cdiagonal=dvector(1,number);


	for(i=1;i<=number;i++)
		X[i]=0;  //set random number vector == 0



//Generate vector of 10 variables ~ N(0,1)

	//fp = fopen("NoiseII.dat","at");

	for(j=1;j<number+1;j++){

		pt1 = ran2(&seed);

		pt2 = ran2(&seed);

		//Z[j]=(sqrt(-2*log(pt1))*cos(2.0*PI*pt2));
		Z[j] = gsl_ran_gaussian(r2,1.0);
		//fprintf(fp,"%d\t%f\n",j,Z[j]);

	}
	//fclose(fp);



//Find the Cholesky Decomposition (C) of Sigma, and Transpose it.



	choldc(Sigma,number,Cdiagonal);  // gives C as lower triangle of Sigma with Cdiagonal with upper triangle == 0 

	

	// Write out C     

	for(j=1;j<number+1;j++){

		for(k=1;k<j;k++){

			C[j][k]=Sigma[j][k];  //take lower triangle of Sigma and make it lower triangle of C

			C[k][j]=0;  	//make upper triangle of C =0

		}

		C[j][j]=Cdiagonal[j];   //diagonal of C

	}





//Multiply C by Vector Z, and add mean to vector X

	for(j=1;j<number+1;j++){

		for(k=1;k<number+1;k++){

			X[j]=X[j]+C[j][k]*Z[k];

		}

		X[j]=X[j]+mean;

	}



	free_dvector(Z,1,number);  

	//free_dmatrix(Sigma,1,number,1,number);  

	free_dmatrix(C,1,number,1,number);

	free_dvector(Cdiagonal,1,number);

	return; 

}

void DistCorrldRand2(float corrln, float mean, float SD, long seed, int MaxRow, int MaxCol, float RowScale, float ColScale, float DistParam, double **Sigma, double *X,double *GaussRandNums){


// Choose Z1,Z2... from N(0,1), Convert to Correl. Normal. random numbers  X1,X2,... 

// by X = C Z  where C is cholesky decomp of the Covariance matrix of X


	/*
	double corrln = 0.999;
	double mean=0, var=0.25, covar;
	long seed = -1;
	*/

	FILE *fp;

	int i,j,k,l;
	int Pop1, Pop2;
	int number = MaxRow*MaxCol;

	//int number=100; // number of random numbers to generate

	double pt1,pt2,*Z, **C,*Cdiagonal;
	float covar, var;
	float Dist;
	long junk;


      
	var = SD*SD;


	covar = corrln*var;

	//X=dvector(1,number);  // Vector of 10 correlated normal random numbers

	Z=dvector(1,number);   // Vector of 10 normal random numbers

	//Sigma=dmatrix(1,number,1,number);  //Matrix of variances and covariances 

	C=dmatrix(1,number,1,number);

	Cdiagonal=dvector(1,number);


	for(i=1;i<=number;i++)
		X[i]=0;  //set random number vector == 0



//Generate vector of 10 variables ~ N(0,1)

	//fp = fopen("NoiseII.dat","at");

	for(j=1;j<number+1;j++){

		pt1 = ran2(&seed);

		pt2 = ran2(&seed);

		//Z[j]=(sqrt(-2*log(pt1))*cos(2.0*PI*pt2));
		Z[j] = GaussRandNums[j];

		//fprintf(fp,"%d\t%f\n",j,Z[j]);

	}
	//fclose(fp);



//Find the Cholesky Decomposition (C) of Sigma, and Transpose it.



	choldc(Sigma,number,Cdiagonal);  // gives C as lower triangle of Sigma with Cdiagonal with upper triangle == 0 

	

	// Write out C     

	for(j=1;j<number+1;j++){

		for(k=1;k<j;k++){

			C[j][k]=Sigma[j][k];  //take lower triangle of Sigma and make it lower triangle of C

			C[k][j]=0;  	//make upper triangle of C =0

		}

		C[j][j]=Cdiagonal[j];   //diagonal of C

	}





//Multiply C by Vector Z, and add mean to vector X

	for(j=1;j<number+1;j++){

		for(k=1;k<number+1;k++){

			X[j]=X[j]+C[j][k]*Z[k];

		}

		X[j]=X[j]+mean;

	}

	if(SD<=0)
	  for(j=1;j<number+1;j++)
	    X[j] = 0.0; 



	free_dvector(Z,1,number);  

	//free_dmatrix(Sigma,1,number,1,number);  

	free_dmatrix(C,1,number,1,number);

	free_dvector(Cdiagonal,1,number);


	return; 

}

void CorrldRand(float corrln, float mean, float SD, long seed, int number, double *X){


// Choose Z1,Z2... from N(0,1), Convert to Correl. Normal. random numbers  X1,X2,... 

// by X = C Z  where C is cholesky decomp of the Covariance matrix of X


	/*
	double corrln = 0.999;
	double mean=0, var=0.25, covar;
	long seed = -1;
	*/

	FILE *fp;

	int j,k;

	//int number=100; // number of random numbers to generate

	double pt1,pt2,*Z, **Sigma, **C,*Cdiagonal;
	float covar, var;
	var = SD*SD;

	
	covar = corrln*var;

	//X=dvector(1,number);  // Vector of 10 correlated normal random numbers

	Z=dvector(1,number);   // Vector of 10 normal random numbers

	Sigma=dmatrix(1,number,1,number);  //Matrix of variances and covariances 

	C=dmatrix(1,number,1,number);

	Cdiagonal=dvector(1,number);


//Define the Covariance matrix   *** Covar[j][k] <=  sqrt(Var[j]) sqrt(Var [k])

// **Must be Symmetric, Positive-Definite Matrix

	for(j=1;j<number+1;j++){

		for(k=1;k<number+1;k++){

			Sigma[j][k]=covar;  //Covariance between random variables

			if(j==k){ Sigma[j][j]=var;}  //Variance within a random variable

		}



		X[j]=0;  //set random number vector == 0

	}





//Generate vector of 10 variables ~ N(0,1)

	for(j=1;j<number+1;j++){

		pt1 = ran2(&seed);

		pt2 = ran2(&seed);

		//Z[j]=(sqrt(-2*log(pt1))*cos(2.0*PI*pt2));
		Z[j] = gsl_ran_gaussian(r2,1.0);

	}





//Find the Cholesky Decomposition (C) of Sigma, and Transpose it.



	choldc(Sigma,number,Cdiagonal);  // gives C as lower triangle of Sigma with Cdiagonal with upper triangle == 0 

	

	// Write out C     

	for(j=1;j<number+1;j++){

		for(k=1;k<j;k++){

			C[j][k]=Sigma[j][k];  //take lower triangle of Sigma and make it lower triangle of C

			C[k][j]=0;  	//make upper triangle of C =0

		}

		C[j][j]=Cdiagonal[j];   //diagonal of C

	}





//Multiply C by Vector Z, and add mean to vector X

	for(j=1;j<number+1;j++){

		for(k=1;k<number+1;k++){

			X[j]=X[j]+C[j][k]*Z[k];

		}

		X[j]=X[j]+mean;

	}



	free_dvector(Z,1,number);  

	free_dmatrix(Sigma,1,number,1,number);  

	free_dmatrix(C,1,number,1,number);

	free_dvector(Cdiagonal,1,number);

	return; 

}

