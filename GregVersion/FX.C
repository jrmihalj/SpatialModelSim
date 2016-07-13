double fx(double k, double oldN, double oldZ, double I)
{
	double dummy;
	
	dummy = (1.0/k)*(oldN*I + oldZ);
	return(1.0 - I - pow((1+dummy),-k));


}
