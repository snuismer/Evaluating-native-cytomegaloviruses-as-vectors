#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>
#include <ctime>
#include <random>  

using namespace std;

time_t seconds;

int replicate,SimSets,RunType;//Variables to set up a multi-run framework
int Step,MaxSteps,Burn;//Variables for MCMC
int test,ani;//Debugging variables only
int Thinner,HowThin; //Thinning for MCMC chain output
int Rows;//The number of rows in the data frame where each row is a unique animal measurement. Calculated upon data ingestion 
double Shed_Data[5000][2];//A container for the data where columns are age and infection status with an individual animal in each row
int Accepts[20],Declines[20];//Containers for acceptance and decline rates by chain

//Define parameters
int Dimension;//Defines the dimension of the model with respect to parameters being estimated 
int NumChain,Chain;//How many chains to run
double VCV[500][500];//The Variance-Covariance matrix for the parameters
double bM,bS,nM,nS,R0M,R0S,R0Lo,R0Hi,omegaM,omegaS,rhoM,rhoS,sigmaAlpha,sigmaBeta,sigLo,sigHi;//The parameters that define the prior distribution
double bC[20],bP[20],muC[20],muP[20],nC[20],nP[20],R0C[20],R0P[20],omegaC[20],omegaP[20],rhoC[20],rhoP[20],sigmaC[20],sigmaP[20];//The current and proposed parameter values 
double PriorCurrent[20],PriorProposed[20];//The prior probability of the current position and proposed position with respect to chain
double LikeCurrent[20],LikeProposed[20];//The likelihood of the data given the current position and a proposed position with respect to chain
double ABar,ACount,IBar,ICount;//Variables for calculating average age and prevalence in the sample
double GRStatR0;//Convergence statistic

ifstream in_Shed;//This is the file stream used to read in the age and infection data
ofstream out_Chains;//This is the file stream used to output the MCMC chains
ofstream out_Converge;//This is the file stream used to output the convergence status of each run (used for simulation study)
ofstream out_Data;//This is the file stream used to output the data file to check that it was ingested correctly


random_device rd;   // non-deterministic generator  
mt19937 gen(rd());  // seed mersenne twister.  

//Define some useful random number generators
double Big_Rand(double MIN, double MAX)
{
	uniform_real_distribution<> dist(MIN, MAX); // distribute results between Min and Max
	return dist(gen);
}

int Rand_Int(int MIN, int MAX)
{
	uniform_int_distribution<> dist(MIN, MAX); // distribute results between Min and Max
	return dist(gen);
}

int Bin_Rand(int Tries, double Srate)
{
	binomial_distribution<> distr(Tries,Srate); // distribute results between 1 and 6 inclusive.  
	return distr(gen);
}

double Norm_Rand(double Mu, double SD)
{ 
	normal_distribution<> distr(Mu, SD); // 
	return distr(gen);
}

double LogNorm_Rand(double p, double q)
{ 
	lognormal_distribution<> distr(p, q); // p is mean and q is sd
	return distr(gen);
}

double Gam_Rand(double p, double q)
{
	gamma_distribution<> distr(p, q); // p is Alpha (shape) and q is Beta (rate)
	return distr(gen);
}

double Exp_Rand(double p)
{
	exponential_distribution<> distr(p); // p is rate
	return distr(gen);
}
//Done defining random number generators

//Define some useful functions
//This function calculates the log-likelihood of observing infection status X given age a as a function of model parameters
//The expressions in this function were derived in the accompanying Mathematica notebook
double LogLike(double b, double R0, double omega, double rho, double sigma, double a, int X)
{
	double L,k1,k2,Force,LV1,LV2,Guts,GutNum;

	//First, check that the value of R_0>1. Values less than 1 are impossible given the model assumptions (steady state viral presence). Practically speaking, we need to treat values of R0<=1 manually or our expressions blow up
	if(R0<=1)
	{
		L=-pow(10,100.0);//Approximate the log of zero with an arbitrarily small number
	}

	if(R0>1)
	{
		//Next, calculate the force of infection at steady state
		Force=(b*(-1 + R0)*(b + omega + rho - (b + rho)*sigma))/(b + omega + rho);
		
		if(X==0)//if the animal is CMV negative, use this
		{
			L=log(1 - (exp(a*Force)*rho*(-Force + omega + rho) + 
		((b + Force)*(Force - rho)*(omega + rho)*(-omega + (b + rho)*(-1 + sigma)))/
			(-((b + Force)*omega) + (b + rho)*(-Force + b*(-1 + sigma))) - 
		(exp(a*(Force - omega - rho))*Force*omega*(b + omega + rho)*(b + Force - (b + rho)*sigma))/
			((b + Force)*omega + (b + rho)*(b + Force - b*sigma)))/(exp(a*Force)*(omega + rho)*(-Force + omega + rho)));
		}
		if(X==1)//if the animal is CMV positive, use this
		{
			L=log((exp(a*Force)*rho*(-Force + omega + rho) + 
		((b + Force)*(Force - rho)*(omega + rho)*(-omega + (b + rho)*(-1 + sigma)))/
		(-((b + Force)*omega) + (b + rho)*(-Force + b*(-1 + sigma))) - 
		(exp(a*(Force - omega - rho))*Force*omega*(b + omega + rho)*(b + Force - (b + rho)*sigma))/
		((b + Force)*omega + (b + rho)*(b + Force - b*sigma)))/(exp(a*Force)*(omega + rho)*(-Force + omega + rho)));
		}
	}
	return L;
}

//The following functions are used to calculate prior probabilities
double LogPriorProbLN(double Z, double X, double Y)//The log probability of drawing Z from a lognormal PDF given X is the MODE and Y is the standard deviation
{
	double PofZ;
	PofZ=-(pow(-pow(Y,2.0)-log(X)+log(Z),2.0))/(2*pow(Y,2))-log(sqrt(2*M_PI)*Y*Z);
	return (PofZ);
}

double LogPriorProbUni(double Z, double X, double Y)//The log probability of drawing Z from a uniform distribution between X and Y
{
	double PofZ;
	if(Z>=X&&Z<=Y)
	{
		PofZ=log(1.0/(Y-X));
	}
	else
	{
		PofZ=-pow(10,100.0);//approximate the log of zero with an arbitrarily small number
	}
	return (PofZ);
}
double LogPriorProbEx(double Z, double X)//The log probability of drawing Z from an exponential PDF given rate parameter X
{
	double PofZ;
	PofZ=log(X)-Z*X;
	return (PofZ);
}
double LogPriorProbBe(double Z, double X, double Y)//The log probability of drawing Z from a Beta PDF with parameters X and Y
{
	double PofZ;
	if(Z>0&&Z<1)
	{
		PofZ=(-1 +Y)*log(1 - Z) + (-1 + X)*log(Z) - log(beta(X,Y));
	}
	else
	{
		PofZ=-pow(10,100.0);//approximate the log of zero with an arbitrarily small number
	}
	return (PofZ);
}

//Define void functions. Names should be self-explanatory
void DefineParms();
void IngestData();
void Initialize();
void Propose();
void Likelihood();
void PriorProb();
void UpdateChain();
void OutChain();
void CheckConverge();


int main()
{
	cout<<"Perform analysis of a single data set (1) or multiple indexed sets (2)\n";//Option to analyze a single data set or loop through multiple data sets where each has a numerically indexed file name
	cin>>RunType;
	if(RunType==1)
	{
		SimSets=1; //Execute a single round of inference
	}
	if(RunType==2)
	{
		cout<<"Enter the number of datasets\n";
		cin>>SimSets; //Execute "SimSets" rounds of inference for each individually named and numerically indexed data set
	}
	
	string InputFile;
	string OutputFile;
	if(RunType==1)//If running a single data set, define the filenames here
	{
		//To run the CMV2 data, use this:
		//InputFile="CMV2_Age_Infection.csv";
		//OutputFile="CMV2_Posterior.csv";

		//To run the CMV3 data, use this:
		InputFile="CMV3_Age_Infection.csv";
		OutputFile="CMV3_Posterior.csv";


	}
	if(RunType==2)//If running multiple data sets, define output filename here (the input filenames are set at run by appending a numerical index to a root filename)
	{
		OutputFile="Test_5_13_25.csv";
	}


	out_Chains.open("Chains_"+OutputFile);
	out_Data.open("Data_"+OutputFile);
	out_Converge.open("Convergence_"+OutputFile);

	out_Converge<<"Replicate,"<<"Convergence\n";

	for(replicate=0;replicate<SimSets;replicate++)//For all available data sets, execute inference
	{
		//Open the proper input file for this replicate
		if(RunType==2)//If analyzing simulated data, open the file with name indexed by replicate
		{
			InputFile="SimData"+std::to_string(replicate)+".csv";
		}
		in_Shed.open(InputFile);

		NumChain=5;//Set the number of chains to run
		IngestData(); //Read in the data file
		DefineParms();//Set priors
		Initialize();//Initialize the simulation

		out_Chains<<"Data_Set,Chain,Step,b,R_0,omega,rho,sigma\n";//Set the header for the chain output file
		Burn=300000;//Set the burn in period
		Thinner=0;//Start the thinning counter at 0
		HowThin=50;//Set the thinning parameter
		MaxSteps=1000000;//Set the number of steps over which to run the MCMC

		for(Chain=0;Chain<NumChain;Chain++)//Initialize the acceptance rate monitors
		{
			Accepts[Chain]=0;
			Declines[Chain]=0;
		}

		for(Step=0;Step<MaxSteps;Step++)//Begin MCMC for this data set
		{
			cout<<replicate<<","<<Step<<"\n";//Send aa progress monitor to the screen. Can be commented out to increase speed

			Propose();
			Likelihood();
			PriorProb();
			UpdateChain();
			OutChain();

			if(Step>MaxSteps-100001)//Start checking convergence in the final 100000 steps. This is for use in the simulation study where we want to code each run by convergence status
			{
				CheckConverge();
			}
			Thinner++;
		}//End MCMC for this data set

	in_Shed.close();//Close the shedding data input file each round so the next file can be opened without conflict (when analyzing multiple data sets)
	if(GRStatR0<1.1)
	{
		out_Converge<<replicate<<","<<1<<"\n";
	}
	else
	{
		out_Converge<<replicate<<","<<0<<"\n";
	}

	out_Converge << std::flush; 

	}//End of multi data set analysis loop

	out_Chains.close();//Close the files
	out_Data.close();
	out_Converge.close();
	
}

void DefineParms()//This function defines the prior distribution
{
	//Set the mode and variance for the b prior
	bM=1/ABar;
    bS=0.2;

	//We assume a lognormal prior for R_0. We set the mode to a rough guess based on prevalence.
	R0M=1.0/(1.0-IBar);
	R0S=1.0;

	omegaM=20.0;//We assume an exponential distribution with this rate for omega
    rhoM=20.0;//We assume an exponential distribution with this rate for rho

	//We assume a Beta prior for sigma because it can be only between 0 and 1
	sigmaAlpha=1.0;
	sigmaBeta=10.0;
}

void IngestData()//This function imports the data from a .csv file (the file must have two columns, age and infection status, with each row representing observations/properties of a single animal)
{
	int i, Row, Column;
	char ch;

	Rows=0;//This is the global variable that tracks the total number of rows (animals) in the data set
	i = 0;
	Row = 0;
	Column = 0;
	while ( !in_Shed.eof() ) // keep reading until end-of-file
	{ 
		in_Shed >> Shed_Data[Row][Column]; //Read the next number into the Array Shed_Data from the csv file "in_Shed"
		if(Column!=1)
		{
			in_Shed>>ch;//read and ignore the comma if you are not in the last column
		}

		i++;
		Column = i % 2;//This modulus operator indicates when to switch to a new row and thus represents the number of columns (e.g., 2 indicates 2 columns)
	
		if (i > 0 && Column == 0)
		{
			Row++;//Increment the local counter for the current row being read in from the file
			Rows++;//Increment the global counter for the size of the data set
		}
  	}
   	in_Shed.close();
    cout << "End-of-file reached.." << endl;

	//Calculate average age of animals for use in setting the prior for b
	ABar=0;
	ACount=0;
	for(Row=0;Row<Rows;Row++)
	{
		ABar=ABar+Shed_Data[Row][0];
		ACount++;
	}
	ABar=ABar/(1.0*ACount);

	//Calculate prevalence for use in setting the prior for R0
	IBar=0;
	ICount=0;
	for(Row=0;Row<Rows;Row++)
	{
		IBar=IBar+Shed_Data[Row][1];
		ICount++;
	}
	IBar=IBar/(1.0*ICount);

	//Output the data so that in can be checked for correct ingestion (can be suppressed)
	out_Data<<"Age,CMV\n";
	for(Row=0;Row<Rows;Row++)
	{
        for(i=0;i<2;i++)
        {
			out_Data<<Shed_Data[Row][i]<<","<<Shed_Data[Row][i];
        }
		out_Data<<"\n";
	}
	//Ingested data output complete

	//Set dimensions of the system
	Dimension=5;//The dimension is equal to the number of parameters to be estimated

}
void Initialize()//This function initializes the simulation
{
	int i,j;
	double bMean,R0Mean,omegaMean,rhoMean,sigmaMean;
	double SumLike,MyV,MyA,MyB,MyAge,MyCMV;
	double ThisBit;

	//Initialize the VCV which defines the multivariate normal used to draw proposed moves. Note that we assume no covariance.
	for (i = 0; i < Dimension; i++) 
	{
    	for (j = 0; j < Dimension; j++) 
		{
			if(i==j)//This is a diagonal and thus a variance
			{
				
				if(i==0)//b
				{
					VCV[i][j]=0.000001;
				}
				if(i==1)//R_0
				{
					VCV[i][j]=5.0;
				}
				if(i==2)//omega
				{
					VCV[i][j]=0.001;
				}
				if(i==3)//rho
				{
					VCV[i][j]=0.001;
				}
				if(i==4)//sigma
				{
					VCV[i][j]=0.01;//0.02;
				}
			}
			else
			{
				VCV[i][j]=0;//If it is an off-diagonal it is a covariance and we set these to zero
			}
		}
	}
	for(Chain=0;Chain<NumChain;Chain++)
	{
		do//Draw starting values at random until you get a vector of parameters that yields a feasible likelihood (i.e., not - inf)
		{
			//Draw the starting value of b in an informed way from the prior
			//First we need to move back into the mean space from the mode space for the parameters with lognormal distributions...
			bMean=log(bM)+pow(bS,2.0);//mode to mean
			//Then we can draw from lognormal to get the starting values for birth rate
			bC[Chain]=LogNorm_Rand(bMean, bS);

			//We draw the strating value of R_0 from a uniform on a plausible range
			R0C[Chain]=Big_Rand(2, 10); //Used in test4

			//The starting values of omega and rho are drwn from exponential distributions. 
			omegaC[Chain]=Exp_Rand(omegaM);
			rhoC[Chain]=Exp_Rand(rhoM);

			//Draw sigma from a uniform distribution
			sigmaC[Chain]=Big_Rand(0,0.2);

			//Check the likelihood of the starting point to make sure it is feasible. Draw at random until the starting point yields a Real likelihood
			SumLike=0;
			for(i=0;i<Rows;i++)//for the number of rows in the data frame 
			{
				MyAge=Shed_Data[i][0]; 
				MyCMV=Shed_Data[i][1];

				ThisBit=LogLike(bC[Chain],R0C[Chain],omegaC[Chain],rhoC[Chain],sigmaC[Chain], MyAge, MyCMV);
				SumLike=SumLike+ThisBit;
			}
		}
		while(isinf(SumLike)==true||isnan(SumLike)==true);
	}//Done looping over chains
}

void Propose()//Propose a move from a multivariate normal
{
	int NonPos;//Used to prevent non-positive proposed moves
	int i,j,k;
    double PropVec1[100],PropVec2[100];
	double LT[500][500];//A lower triangular matrix for the Cholesky decomposition
	double sum;

	//The code is set up so that the proposal can be drawn from a multivariate gaussian distribution. Because we assume each dimension is independent, this is not neccessary for this project as described in the accopanying publication
	//Although maintaining the multivariate generality slows the computations, this is not an issue as the procedure is sufficiently rapid even with this extra step
    //We begin by performing the Cholesky decomposition
	for (i = 0; i < Dimension; i++) 
	{
    	for (j = 0; j <= i; j++) 
		{
			LT[i][j]=0;
		}
	}
	//Calculate the LT matrix using the Cholesky Decomposition
	for (i = 0; i < Dimension; i++) 
	{
    	for (j = 0; j <= i; j++) 
		{
			sum = 0;
			for (k = 0; k < j; k++)
			{
				sum += LT[i][k] * LT[j][k];
			}

			if (i == j)
			{
				LT[i][j] = sqrt(VCV[i][i] - sum);
			}
			else
			{
				LT[i][j] = (1.0 / LT[j][j] * (VCV[i][j] - sum));
			}
   		 }
	}

	for(Chain=0;Chain<NumChain;Chain++)
	{
		do//keep drawing if any proposed moves are non-positive
		{
			NonPos=0;
			//First draw a vector of univariate random normals with mean zero and unit deviation. I assume this is a column vector.
			for(i=0;i<Dimension;i++)
			{
				PropVec1[i]=Norm_Rand(0, 1);
				PropVec2[i]=0.0;
			}
			//Then transform these into draws from a multivariate normal using the Cholesky
			for (i = 0; i < Dimension; i++) 
			{
				for (j = 0; j < Dimension; j++) 
				{
					PropVec2[i]=PropVec2[i]+PropVec1[j]*LT[i][j];//Matrix multiplication to yield final vector PropVec2
				}
			}

			//Now, use these deviates to create a new vector of proposed parameters
			bP[Chain]=bC[Chain]+PropVec2[0];
			if(bP[Chain]<=0){NonPos++;}
			R0P[Chain]=R0C[Chain]+PropVec2[1];
			if(R0P[Chain]<=0){NonPos++;}
			omegaP[Chain]=omegaC[Chain]+PropVec2[2];
			if(omegaP[Chain]<=0){NonPos++;}
			rhoP[Chain]=rhoC[Chain]+PropVec2[3];
			if(rhoP[Chain]<=0){NonPos++;}
			sigmaP[Chain]=sigmaC[Chain]+PropVec2[4];
			if(sigmaP[Chain]<=0){NonPos++;}
		}
		while(NonPos>0);//Continue drawing until all proposed parameter values are positive
	}//end loop over chains

}

void Likelihood()//Calculate the likelihood for current and proposed positions
{
	int i,j,k,iMax,jMax,kMax;
	double SumLike;
	double MyAge,MyCMV;

	for(Chain=0;Chain<NumChain;Chain++)
	{
		//First do the likelihood at the current position
		SumLike=0;
		for(i=0;i<Rows;i++)//For each animal in the data set
		{
			MyAge=Shed_Data[i][0]; 
			MyCMV=Shed_Data[i][1];

			SumLike=SumLike+LogLike(bC[Chain],R0C[Chain],omegaC[Chain],rhoC[Chain],sigmaC[Chain],MyAge, MyCMV);
			//cout<<"SumLike= "<<SumLike<<"\n";
		}
		LikeCurrent[Chain]=SumLike;

		if(isnan(LikeCurrent[Chain])==true)//Warn and wait
		{
			cout<<"WARNING: You have a NAN current likelihood\n";
			cin>>test;
		}
		if(isinf(LikeCurrent[Chain])==true)//Warn and wait
		{
			cout<<"WARNING: You now have a INF current likelihood\n";
			cin>>test;
		}
		

		//Then do the likelihood for the proposed position
		SumLike=0;
		for(i=0;i<Rows;i++)//For each animal in the data set
		{
			MyAge=Shed_Data[i][0]; 
			MyCMV=Shed_Data[i][1];

			SumLike=SumLike+LogLike(bP[Chain],R0P[Chain],omegaP[Chain],rhoP[Chain],sigmaP[Chain],MyAge, MyCMV);
		}
		LikeProposed[Chain]=SumLike;
	}

}

void PriorProb()//Calculate prior probability for the current and proposed position
{
	int i;
	double SumProd;

	for(Chain=0;Chain<NumChain;Chain++)
	{
		//First calculate the prior probability for the current position
		SumProd=0;//Start with probability 0
		SumProd=SumProd+LogPriorProbLN(bC[Chain], bM, bS);
		SumProd=SumProd+LogPriorProbLN(R0C[Chain], R0M, R0S);
		SumProd=SumProd+LogPriorProbEx(omegaC[Chain], omegaM);
		SumProd=SumProd+LogPriorProbEx(rhoC[Chain], rhoM);
		SumProd=SumProd+LogPriorProbBe(sigmaC[Chain], sigmaAlpha, sigmaBeta);
		PriorCurrent[Chain]=SumProd;

		//Second calculate the prior probability for the proposed position
		SumProd=0;//Start with probability 0
		SumProd=SumProd+LogPriorProbLN(bP[Chain], bM, bS);
		SumProd=SumProd+LogPriorProbLN(R0P[Chain], R0M, R0S);
		SumProd=SumProd+LogPriorProbEx(omegaP[Chain], omegaM);
		SumProd=SumProd+LogPriorProbEx(rhoP[Chain], rhoM);
		SumProd=SumProd+LogPriorProbBe(sigmaP[Chain], sigmaAlpha, sigmaBeta);		
		PriorProposed[Chain]=SumProd;
	}
}

void UpdateChain()//Decide whether to accept or reject the proposed move
{
	int i,j,k,n,counter,TryAgain;
	double P_Accept,Tester,DifLogLike,DifLogPrior,RatPrior,RatLike,JustRat,SumRat;
	double XSum[500],X2Sum[500][500];//A storage container and summing tools for expected values and expected products

	//For a move to be accepted we require that the liklihood of the proposed move did not evaluate to a NAN in addition to the standard criterion
	for(Chain=0;Chain<NumChain;Chain++)
	{
		//first, check that we have not proposed an infinite or nan likelihood. If we have, preclude acceptance
		TryAgain=0;
		if(isinf(LikeProposed[Chain])==true)
		{
			cout<<"WARNING: Your proposed likelihood is INF\n";
			TryAgain=1;
		}
		
		if(isnan(LikeProposed[Chain])==true)
		{
			cout<<"WARNING: Your proposed likelihood is NAN\n";
			TryAgain=1;
		}

        //Calculate the difference in likelihood and prior probability for current and proposed position and then sum them
		DifLogLike=LikeProposed[Chain]-LikeCurrent[Chain];
		DifLogPrior=PriorProposed[Chain]-PriorCurrent[Chain];
		SumRat=DifLogLike+DifLogPrior;

		Tester=log(Big_Rand(0,1));
		if(Tester<SumRat&&TryAgain==0)//If this is true, accept the proposal
		{
			Accepts[Chain]++;//Increment the rolling average acceptance rate monitor

			//Update values
			bC[Chain]=bP[Chain];
			R0C[Chain]=R0P[Chain];
			omegaC[Chain]=omegaP[Chain];
			rhoC[Chain]=rhoP[Chain];
			sigmaC[Chain]=sigmaP[Chain];
		}
		else//Otherwise, record as a failure by incrementing the declination rate monitor
		{
			Declines[Chain]++;
		}
	}//End chain loop

}

void OutChain()//Output the chains if we are beyond the burn threshold
{

	if(Step==Burn)//Once we hit the burn in threshold, set the thinning counter to zero
	{
		Thinner=0;
	}
	if(Step>Burn&&Thinner==HowThin)//If we are post-burn and at the proper thinning level, output the chain
	{
		for(Chain=0;Chain<NumChain;Chain++)
		{
			out_Chains<<replicate<<","<<Chain<<","<<Step<<","<<bC[Chain]<<","<< R0C[Chain]<<","<< omegaC[Chain]<<","<<rhoC[Chain]<<","<<sigmaC[Chain]<<"\n";
		}
		Thinner=0;
	}
}

void CheckConverge()//Perform convergence checks over the final 100,000 steps to determine convergence status (for use in simulation study)
{
	double L,SumR0[20],SumR02[20],SumER0,SumER02,ER0[20],ER02[20],VR0[20],WCVarR0,BCVarR0;
	double GRStatOmega,SumOmega[20],SumOmega2[20],SumEOmega,SumEOmega2,EOmega[20],EOmega2[20],VOmega[20],WCVarOmega,BCVarOmega;

	if(Step==MaxSteps-100000)//If we are in the final 100,000 steps, check convergence
	{
		for(Chain=0;Chain<NumChain;Chain++)
		{
			SumR0[Chain]=0;
			SumR02[Chain]=0;
			SumOmega[Chain]=0;
			SumOmega2[Chain]=0;
		}
	}
	
	L=1.0*(100000);
	WCVarR0=0;
	WCVarOmega=0;
	for(Chain=0;Chain<NumChain;Chain++)
	{
		SumR0[Chain]=SumR0[Chain]+R0C[Chain];
		SumR02[Chain]=SumR02[Chain]+R0C[Chain]*R0C[Chain];
		SumOmega[Chain]=SumOmega[Chain]+omegaC[Chain];
		SumOmega2[Chain]=SumOmega2[Chain]+omegaC[Chain]*omegaC[Chain];
	
		ER0[Chain]=SumR0[Chain]/(L);
		ER02[Chain]=SumR02[Chain]/(L);
		VR0[Chain]=ER02[Chain]-ER0[Chain]*ER0[Chain];//Calculate variance within the chain
		WCVarR0=WCVarR0+VR0[Chain];//Sum up WI chain variance
		
		EOmega[Chain]=SumOmega[Chain]/(L);
		EOmega2[Chain]=SumOmega2[Chain]/(L);
		VOmega[Chain]=EOmega2[Chain]-EOmega[Chain]*EOmega[Chain];//Calculate variance within the chain
		WCVarOmega=WCVarOmega+VOmega[Chain];//Sum up WI chain variance
	}
	WCVarR0=WCVarR0/(1.0*NumChain);//Expected within chain variance across chains for R_0
	WCVarOmega=WCVarOmega/(1.0*NumChain);//Expected within chain variance across chains for omega

	SumER0=0;
	SumER02=0;
	SumEOmega=0;
	SumEOmega2=0;
	for(Chain=0;Chain<NumChain;Chain++)
	{
		SumER0=SumER0+ER0[Chain];
		SumER02=SumER02+ER0[Chain]*ER0[Chain];
		SumEOmega=SumEOmega+EOmega[Chain];
		SumEOmega2=SumEOmega2+EOmega[Chain]*EOmega[Chain];
	}
	BCVarR0=SumER02/(1.0*(NumChain))-SumER0*SumER0/(1.0*(NumChain)*(NumChain));
	BCVarOmega=SumEOmega2/(1.0*(NumChain))-SumEOmega*SumEOmega/(1.0*(NumChain)*(NumChain));

	GRStatR0=(((L-1.0)/(L))*WCVarR0+(1.0/L)*BCVarR0)/(WCVarR0);
	GRStatOmega=(((L-1.0)/(L))*WCVarOmega+(1.0/L)*BCVarOmega)/(WCVarOmega);
}