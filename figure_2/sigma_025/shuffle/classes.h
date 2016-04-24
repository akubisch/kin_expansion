const int 	XMAX 	= 200;				//max. x-dimension of the world
const int 	YMAX	= 100;				//max. y-dimension of the world
const float VARIANCE	= 0.2;			//variance used for many calculations
const float	MAX_MUT 	= 0.05;			//maximum possible value, that is added by mutation to
//to the dispersal propensity alleles
const short int	MAXF	= 500;
const short int	MAXM	= 500;
const float SD_LD 	= 0.2;          	//standard deviations for gaussan distributed variation

const char  TIME_OUT_INT = 10;

const int   BURN_IN_TIME = 1000;
const int   BURN_IN_REGION = 50;             //how many stripes


class tInd
{
public:
	double  ld[2];						//one locus with 2 dispersal alleles
};

class tSpecies
{
public:
	tSpecies();
	vector<tInd> males,females;			//population
	vector<tInd> malimmis,femimmis;		//immigrants
};

tSpecies::tSpecies() {
	males.clear();
	males.reserve(MAXM);
	females.clear();
	females.reserve(MAXF);
	malimmis.clear();
	malimmis.reserve(MAXM);
	femimmis.clear();
	malimmis.reserve(MAXF);
}

class tParameters
{
public:
	int 	capacity;					//habitat capacity
	double	lambda_0;					//per capita growth rate
	double	disp_mort;					//dispersal mortality
	double	sigma;						//environmental fluctuations
	double	ext_prob;					//extinction rate
	double	mut_prob;					//mutation probability
	int		tmax;						//max. no. of generations to be simulated

	int replications;					//number of replicates per experiment

	bool    burnin;                     //is currently the burn-in running?
};

class tAnalysis
{
public:
	tAnalysis();
	vector<unsigned short int> 			times;					//stores all generation times
	vector<unsigned long> 	N;					//saves metapopulation sizes of species 1
	vector<double>			emi;					//emigration rate of species 2
	vector<double>			emi_cor;
	vector<double>			emi_mar;
	unsigned long			sum_ind_b4;	//individuals, that may disperse
	unsigned long			sum_disp;	//individuals, that have dispersed
	unsigned long     sum_ind_b4_cor;
	unsigned long     sum_ind_b4_mar;
	unsigned long     sum_disp_cor;
	unsigned long     sum_disp_mar;
	vector<double>			real_border;			//realized range border
	vector<double>			real_border_50;			//realized range border
};

tAnalysis::tAnalysis() {
	times.clear();
	N.clear();
	emi.clear();
	emi_cor.clear();
	emi_mar.clear();
	real_border.clear();
	real_border_50.clear();
}
