/*
*
*		     Simulations of Insects in Spatially Changing Environments
*      >> Multi-generation & kin selection during range expansions <<
*	  			  -----------------------------
*		      			Alexander Kubisch
*				   Evolutionary Ecology Group
*			  	    University of Wuerzburg
*
*/

#include <iostream>
#include <cstdlib>							//standard C library
#include <ctime>							//access system time library
#include <fstream>							//file streaming library
#include <string>							//string library included
#include <sstream>							//string streaming for reading numbers from
//infiles
#include <vector>
#include <cmath>							//standard math library

#include <gsl/gsl_rng.h>					//gsl random number generator
#include <gsl/gsl_randist.h>				//gsl random distributions
#include <gsl/gsl_statistics.h>				//gsl some statistical methods
#include <gsl/gsl_statistics_int.h>			//gsl some integer statistical methods
#include <gsl/gsl_sort_int.h>				//gsl sort integer arrays
#include <gsl/gsl_math.h>					//provides additional standard math functions
#include <gsl/gsl_histogram.h>                          //provides histogram functions needed for relatedness

using namespace std;

#include "procedures.h"						//procedure simplifications
#include "classes.h"						//class definitions

//Variables

vector<vector<tSpecies> > 		world;				//simulated landscape
tParameters	par;							//parameters
tAnalysis 	ana;
double maxtemp=0;
double mintemp=1000;
double initial_range;

int			real_border,pred_border;		//realized and predicted range border

fstream 	parinfile("data/para.in", ios::in);		//parameter infile
fstream		outfile_noshuf("data/output_noshuffle.txt", ios::out);	//result outfile without shuffling
fstream		outfile_shuf("data/output_shuffle.txt", ios::out);	//result outfile with shuffling

///////////////////////////////////////////////////////////////////////////////////////
//---------------------------------SET PARAMETERS------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

void set_parameters() {					//read in standard simulation parameters

  string buffer;
  istringstream is;

  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.capacity;					//habitat capacity
  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.lambda_0;					//per capita growth rate
  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.sigma;					//environmental stochasticity
  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.disp_mort;				//dispersal mortality
  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.ext_prob;					//extinction probability
  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.mut_prob;					//mutation probability
  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.tmax;						//simulated time
  getline(parinfile,buffer); is.clear(); is.str(buffer);
  is >> par.replications;				//no of replications
}

///////////////////////////////////////////////////////////////////////////////////////
//-----------------------------INITIALIZATION FUNCTIONS------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

void init_individual(tInd *newind) {
  newind->ld[0] = ran();
  newind->ld[1] = ran();
}

void init_world() {

  tInd newind;

  world.resize(XMAX);
  for (int x=0; x < XMAX; x++) {
    world[x].resize(YMAX);
  }

  for (int x = 0; x < XMAX; x++)
  for (int y = 0; y < YMAX; y++) {

    world[x][y].females.clear();
    world[x][y].males.clear();

    if (x<BURN_IN_REGION) {
      for (int i = 0; i < (par.capacity/2); i++) {
        init_individual(&newind);
        world[x][y].females.push_back(newind);
      }
      for (int i = 0; i < (par.capacity/2); i++) {
        init_individual(&newind);
        world[x][y].males.push_back(newind);
      }
    }

  }

}

///////////////////////////////////////////////////////////////////////////////////////
//-------------------------------SIMULATION FUNCTIONS--------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
//-Dispersal---------------------------------------------------------------------------

double disp_prob(tInd ind) {
  double trait = mean(ind.ld,2);
  return(trait);
}

void find_target_patch(int xs, int ys, int *xt, int *yt) {
  double dist;				//drawn dispersal distance
  double phi;

  dist = 1;						//gives the possibility to implement a kernel

  phi = 2 * M_PI * ran();			//direction is uniformly drawn (0°-360°)

  double xstart = double(xs) + ran() - 0.5;		//simulating area to area dispersal
  double ystart = double(ys) + ran() - 0.5;		//(avoid artefacts by drawing starting position
    //for dispersal randomly inside the patch)

    int xtt = rnd(xstart + dist*cos(phi));
    int ytt = rnd(ystart + dist*sin(phi));


    if (par.burnin) {
      if (xtt>=BURN_IN_REGION) {xtt = 0;}
      if (xtt<0) {xtt = BURN_IN_REGION-1;}
    }

    *xt = xtt;
    *yt = ytt;
  }

  void dispersal(int x, int y, int xt, int yt, unsigned int i, int sex) {

    if (yt<0) {yt = YMAX;}; if (yt>=YMAX) {yt = 0;};	//building a tube

    if (sex==0) {	//female dispersal
      if ((ran()>par.disp_mort)&&	//if survival & target inside the world
      (xt>=0)&&(xt<XMAX)) {					//a tube with absorbing x-conditions
        world[xt][yt].femimmis.push_back(world[x][y].females.at(i));	//immigration
      }
    }
    if (sex==1) {	//male dispersal
      if ((ran()>par.disp_mort)&&
      (xt>=0)&&(xt<XMAX)) {
        world[xt][yt].malimmis.push_back(world[x][y].males.at(i));
      }
    }
  }

  void dispersal_loop(int t) {
    int xt,yt;

    ana.sum_ind_b4 = 0;
    ana.sum_disp = 0;
    ana.sum_ind_b4_cor = 0;
    ana.sum_ind_b4_mar = 0;
    ana.sum_disp_cor = 0;
    ana.sum_disp_mar = 0;

    int range = 0;
    int core = 24;
    if (t>TIME_OUT_INT) {range = round(ana.real_border_50.back())-2;}

    for (int x = 0; x < XMAX; x++)
    for (int y = 0; y < YMAX; y++) {

      int fpop = world[x][y].females.size();
      int mpop = world[x][y].males.size();

      if ((fpop+mpop)>0) {	//dispersal only, when there are beetles in the patch

        for (unsigned int f = 0; f < world[x][y].females.size(); f++) {
          double d = disp_prob(world[x][y].females.at(f));
          ana.sum_ind_b4++;
          if (x==core) {ana.sum_ind_b4_cor++;}
          if (x==range){ana.sum_ind_b4_mar++;}

          if (ran()<d) {
            find_target_patch(x,y,&xt,&yt);
            dispersal(x,y,xt,yt,f,0);		//dispersal
            ana.sum_disp++;
            if (x==core) {ana.sum_disp_cor++;}
            if (x==range){ana.sum_disp_mar++;}
          } else {
            world[x][y].femimmis.push_back(world[x][y].females.at(f));
          }

        }
        for (unsigned int m = 0; m < world[x][y].males.size(); m++) {
          double d = disp_prob(world[x][y].males.at(m));
          ana.sum_ind_b4++;
          if (x==core) {ana.sum_ind_b4_cor++;}
          if (x==range){ana.sum_ind_b4_mar++;}

          if (ran()<d) {
            find_target_patch(x,y,&xt,&yt);
            dispersal(x,y,xt,yt,m,1);
            ana.sum_disp++;
            if (x==core) {ana.sum_disp_cor++;}
            if (x==range){ana.sum_disp_mar++;}
          } else {
            world[x][y].malimmis.push_back(world[x][y].males.at(m));
          }
        }

      }

    }

    for (int x = 0; x < XMAX; x++)
    for (int y = 0; y < YMAX; y++) {
      world[x][y].females 	= world[x][y].femimmis;
      world[x][y].males 	= world[x][y].malimmis;
      world[x][y].femimmis.clear();
      world[x][y].malimmis.clear();
    }

  }

  //_____________________________________________________________________________________
  //-Reproduction------------------------------------------------------------------------

  void genetics(tInd *newborn, tInd mother, tInd father) {

    if (ran()>0.5) {newborn->ld[0] = mother.ld[0];} else {newborn->ld[0] = mother.ld[1];}
    if (ran()>0.5) {newborn->ld[1] = father.ld[0];} else {newborn->ld[1] = father.ld[1];}

    if (ran()<par.mut_prob) {
      newborn->ld[0] += gauss(VARIANCE);
    }
    if (ran()<par.mut_prob) {
      newborn->ld[1] += gauss(VARIANCE);
    }
  }

  void reproduction(int x, int y, double lambda, double survival) {

    vector<tInd> girls,boys;
    tInd newborn;

    girls.clear(); boys.clear();

    for (unsigned int f = 0; f < world[x][y].females.size(); f++) {
      int kids = poisson(2.0*lambda*survival);	//children numbers are poisson distributed
      int partner = floor(ran() * world[x][y].males.size());	//partner is searched at random

      for (int child = 0; child < kids; child++) {
        float sex = ran();	//gender is randomly drawn

        genetics(&newborn,
          world[x][y].females.at(f),
          world[x][y].males.at(partner));
          double surv_prob = 2;//adaptation(x,y,s,newborn);

          if (ran()<surv_prob) {
            if (sex>=0.5) {girls.push_back(newborn);} else {boys.push_back(newborn);}
          }

        }
      }

      world[x][y].females = girls;
      world[x][y].males = boys;
    }

    void reproduction_loop() {

      for (int x = 0; x < XMAX; x++)
      for (int y = 0; y < YMAX; y++) {

        double lambda = lognorm(par.lambda_0,par.sigma);

        unsigned long sum_pop;
        sum_pop = 0;
        sum_pop += world[x][y].females.size()+
        world[x][y].males.size();

        double a = (par.lambda_0-1)/double(par.capacity);
        double survival = 1/(1+a*double(sum_pop));	//assuming logistic growth

        if ((world[x][y].females.size()>0)&&(world[x][y].males.size()>0)) {
          reproduction(x,y,lambda,survival);
        } else {
          world[x][y].females.clear();
          world[x][y].males.clear();
        }

      }
    }

    //_____________________________________________________________________________________
    //-Extinction--------------------------------------------------------------------------

    void extinction() {
      for (int x = 0; x < XMAX; x++)
      for (int y = 0; y < YMAX; y++) {
        if (ran()<par.ext_prob) {
          world[x][y].females.clear();
          world[x][y].males.clear();
        }
      }
    }

    //_____________________________________________________________________________________
    //-Shuffling of Individuals within one stripe------------------------------------------

    void shuffle() {
      for (int x=0; x < XMAX; x++){
        vector<int> pops_f,pops_m;
        vector<int> no_f,no_m;
        vector<tInd> allfemales,allfemales2;
        vector<tInd> allmales,allmales2;
        pops_f.clear(); pops_m.clear();
        no_f.clear(); no_m.clear();
        allfemales.clear(); allmales.clear();
        allfemales2.clear(); allmales2.clear();
        int cnt_f = 0; int cnt_m = 0;

        for (int y=0; y < YMAX; y++){
          pops_f.push_back(world[x][y].females.size());
          pops_m.push_back(world[x][y].males.size());
          for (unsigned int f = 0; f < world[x][y].females.size(); f++){
            allfemales.push_back(world[x][y].females.at(f));
            no_f.push_back(cnt_f);
            cnt_f++;
          }
          for (unsigned int m = 0; m < world[x][y].males.size(); m++){
            allmales.push_back(world[x][y].males.at(m));
            no_m.push_back(cnt_m);
            cnt_m++;
          }
        }

        int cntf[cnt_f];
        int cntm[cnt_m];
        for (int f=0; f<cnt_f; f++){
          cntf[f] = no_f.at(f);
        }
        for (int m=0; m<cnt_m; m++){
          cntm[m] = no_m.at(m);
        }
        if (cnt_f>0) {shuffle_int(cntf,cnt_f);}
        if (cnt_m>0) {shuffle_int(cntm,cnt_m);}

        cnt_f=0; cnt_m=0;
        for (int y=0; y<YMAX; y++){
          world[x][y].females.clear();
          for (int f=0; f<pops_f.at(y); f++) {
            int current = cntf[cnt_f];
            world[x][y].females.push_back(allfemales.at(current));
            cnt_f++;
          }
          world[x][y].males.clear();
          for (int m=0; m<pops_m.at(y); m++) {
            int current = cntm[cnt_m];
            world[x][y].males.push_back(allmales.at(current));
            cnt_m++;
          }
        }
      }

    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //------------------------------ANALYZING FUNCTIONS----------------------------------//
    ///////////////////////////////////////////////////////////////////////////////////////

    //_____________________________________________________________________________________
    //-Analysis----------------------------------------------------------------------------

    void analyze_range_border(double *r50) {

      bool found = false;
      int r50_2;
      double occ[XMAX];
      int occcnt = 0;

      for (int x=0; x < XMAX; x++) {
        occcnt = 0;
        for (int y=0; y < YMAX; y++) {
          if ((world[x][y].females.size()>0)&
          (world[x][y].males.size()>0)) {
            occcnt++;
          }
        }
        occ[x] = double(occcnt)/double(YMAX);
      }

      bool there_is_a_range = false;
      for (int x=0; x < XMAX; x++) {
        if (occ[x]>0.05) {there_is_a_range=true;}
        if ((occ[x-1]>0.05)&(occ[x]<0.05)) {
          r50_2 = x;
          found = true;
        }
        if (found) {break;}
      }

      if ((!found)&(there_is_a_range)) {r50_2 = ana.real_border_50.back();}
      if ((!found)&(!there_is_a_range)) {r50_2 = 0;}
      *r50 = r50_2;
    }

    //_____________________________________________________________________________________
    //-Analyze-----------------------------------------------------------------------------

    void analyze(int t) {
      ana.times.push_back(t);

      unsigned long s,sn;
      s = 0;
      sn = 0;

      int cnt = 0;

      for (int y = 0; y < YMAX; y++)
      for (int x = 0; x < XMAX; x++) {
        s = world[x][y].females.size()+world[x][y].males.size();
        sn += s;

        cnt++;
      }
      ana.N.push_back(sn);

      if (ana.sum_ind_b4>0) {ana.emi.push_back(double(ana.sum_disp)/double(ana.sum_ind_b4));} else {ana.emi.push_back(0);}

      if (ana.sum_ind_b4_cor>0) {ana.emi_cor.push_back(double(ana.sum_disp_cor)/double(ana.sum_ind_b4_cor));} else {ana.emi_cor.push_back(0);}
      if (ana.sum_ind_b4_mar>0) {ana.emi_mar.push_back(double(ana.sum_disp_mar)/double(ana.sum_ind_b4_mar));} else {ana.emi_mar.push_back(0);}

      double r50;	//range margin location and width
      if (ana.N.back()>0) {
        analyze_range_border(&r50);

        ana.real_border_50.push_back(r50);

      } else {
        double p150 = ana.real_border_50.back();
        ana.real_border_50.push_back(p150);
      }

    }

    //_____________________________________________________________________________________
    //-Write out Analysis Results

    void write_results(int r, int t, bool shuf) {
      if (shuf) {
        outfile_shuf << r << "  " << ana.times.at(t)+1 << "  " << "	" << ana.real_border_50.at(t) << "  " << ana.N.at(t) << "  " << ana.emi.at(t) << "  " << ana.emi_cor.at(t) << "  " << ana.emi_mar.at(t) << endl;
      } else {
        outfile_noshuf << r << "  " << ana.times.at(t)+1 << "  " << "	" << ana.real_border_50.at(t) << "  " << ana.N.at(t) << "  " << ana.emi.at(t) << "  " << ana.emi_cor.at(t) << "  " << ana.emi_mar.at(t) << endl;
      }
    }

    //_____________________________________________________________________________________
    //-Clear Analysis----------------------------------------------------------------------

    void clear_analysis() {
      ana.emi.clear();
      ana.emi_cor.clear();
      ana.emi_mar.clear();
      ana.real_border.clear();
      ana.real_border_50.clear();
      ana.N.clear();
      ana.times.clear();
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //----------------------------------MAIN FUNCTIONS-----------------------------------//
    ///////////////////////////////////////////////////////////////////////////////////////

    //_____________________________________________________________________________________
    //-Initialize--------------------------------------------------------------------------

    void initialize() {

      init_world();

      clear_analysis();

    }


    void simulate(bool shuf, int r) {

      // burn in

      for (int t=0; t < BURN_IN_TIME; t++) {

        cout << "burn-in t: " << t << endl;

        par.burnin = true;

        dispersal_loop(t);
        reproduction_loop();
        extinction();

        analyze(t);

        if (((t+1) % TIME_OUT_INT) == 0) {
          write_results(r+1,t,shuf);
        };

      }

      for (int t = BURN_IN_TIME; t < (par.tmax+BURN_IN_TIME); t++) {

        cout << "t: " << t << endl;

        par.burnin = false;

        dispersal_loop(t);
        reproduction_loop();
        extinction();
        if (shuf) {shuffle();}

        analyze(t);

        if (((t+1) % TIME_OUT_INT) == 0) {
          write_results(r+1,t,shuf);
        };

      }
    }

    int main() {

      cout.setf(ios_base::scientific);		//exponential data output

      specify_rng(time(NULL));			//random randomization?

      outfile_noshuf << "repli	t	real_b	pop	emi emi_cor emi_mar" << endl;
      outfile_shuf << "repli	t	real_b	pop	emi emi_cor emi_mar" << endl;


      cout << "simulating control experiment" << endl;

      set_parameters();
      for (int r = 0; r < par.replications; r++) {
        cout << "replication: " << r+1 << endl;
        initialize();
        simulate(true,r);
      }

      outfile_noshuf.close();
      outfile_shuf.close();

      return 0;
    }
