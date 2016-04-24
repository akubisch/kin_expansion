const gsl_rng *gBaseRand;

//________________________________________________________________________________________
//------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = random();
	gsl_rng_set(gBaseRand, r);
}

//________________________________________________________________________________________
//-------------------------------------------------------------------------Simplifications

//------------------------------------------------------------------------Shuffle an Array

void shuffle_int(int array[], int size)
{
        gsl_ran_shuffle(gBaseRand,array,size,sizeof(int));
}

//-------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//------------------------------------------------------------------Simplify Cauchy Randoms

double cauchy(double d)
{
	return gsl_ran_cauchy(gBaseRand,d);
}

//-----------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}

//---------------------------------------------------------------Simplify Lognormal Random

double lognorm(double zeta, double sigma)
{
	double var;											//variance of resulting randoms
	double mu,s_d;										//mean and sd of lognormal distr.
														//to be calculated by mean and
														//sigma of resulting randoms
	var = sigma*sigma;
	s_d = sqrt(log((var/(2*zeta))+1));
	mu = log(zeta)-0.5*(s_d*s_d);
	return gsl_ran_lognormal(gBaseRand,mu,s_d);
}

//---------------------------------------------------------------Simplify mean calculation

double mean(double data[], int n)
{
	return gsl_stats_mean(data,1,n);
}

//---------------------------------------------------------------Simplify mean calculation

double mean(float data[], int n)
{
    double dat[n];
    for(int i=0; i<n; i++) {dat[i] = double(data[i]);}
	return gsl_stats_mean(dat,1,n);
}

//-----------------------------------------------------Simplify integer median calculation

double median(vector<int> vec, int n)
{
	int data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort_int(data,1,n);
	double medi = gsl_stats_int_median_from_sorted_data(data,1,n);
	return(medi);
}

//-----------------------------------------------------Simplify vector mean calculation

double mean(vector<double> vec, int n)
{
        double data[n];
        for(int i = 0; i<n; i++) {
                data[i] = vec.at(i);
        }
        double m = gsl_stats_mean(data,1,n);
        return(m);
}

//---------------------------------------------------Simplify integer quantile calculation

double q005(vector<int> vec, int n)
{
	int data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort_int(data,1,n);
	double q = gsl_stats_int_quantile_from_sorted_data(data,1,n,0.05);
	return(q);
}

double q095(vector<int> vec, int n)
{
	int data[n];
	for(int i = 0; i<n; i++) {
		data[i] = vec.at(i);
	}
	gsl_sort_int(data,1,n);
	double q = gsl_stats_int_quantile_from_sorted_data(data,1,n,0.95);
	return(q);
}

//----------------------------------------------------Simplify minimum and maximum finding

double min_int(vector<int> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = double(vec.at(i));
	}
	return(round(gsl_stats_min(data,1,n)));
}

double max_int(vector<int> vec, int n)
{
	double data[n];
	for(int i = 0; i<n; i++) {
		data[i] = double(vec.at(i));
	}
	return(round(gsl_stats_max(data,1,n)));
}

//----------------------------------------------------------------------------Simplify Sum

int sum_int(int data[], int n)
{
    int s=0;
    for(int i=0; i<n; i++){
        s += data[i];
    }
    return(s);
}

//-------------------------------------------------------------Simplify standard deviation

double sd(double data[], int n)
{
	return gsl_stats_sd(data,1,n);
}

//---------------------------------------------------------------------Round double to int

int rnd(double a)
{
	return round(a);
}

//----------------------------------------------------------------Test equality of doubles

/*
This function compares two double numbers by using their internal representation as integer
numbers, lexicographically ordered to avoid problems with the 0-limit. "maxUlps" determines
the number of doubles that is tolerated and thus being assumed to be equal to each other.
Is hard to explain in short, have a look at the following website for further explanation:
http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
*/

bool equal(double A, double B, int maxUlps)
{
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.

    int aInt = *(int*)&A;

    // Make aInt lexicographically ordered as a twos-complement int
    if (aInt < 0)
        aInt = 0x80000000 - aInt;

    // Make bInt lexicographically ordered as a twos-complement int
    int bInt = *(int*)&B;

    if (bInt < 0)
        bInt = 0x80000000 - bInt;

    int intDiff = abs(aInt - bInt);
    if (intDiff <= maxUlps)
        return true;
    return false;
}
