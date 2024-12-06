#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <iterator>
#include <algorithm>
#include "cminpack.h"
#include <stdlib.h>
#include <iomanip>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
double alpha( double G, const double *pars );
double zbar( double en );
double espec( double e, double a, double b, double emax );
int Chi2( void *p, int m, int n, const double *pars, double *fvec, int iflag);

int nDataPoints;
int num_pars = -1;
int grad_order = -1;
int bootstrap_opt = -1;
int rand_offset = -1;
struct DataPoint{
	double G;
	double AlphaTop;
	double AlphaErrTop;
	double AlphaBot;
	double AlphaErrBot;
};
std::vector<DataPoint> g10_data;
std::vector<DataPoint> bootstrap_g10_data;


// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [bootstrap opt] \n";
		cerr << "**********************************\n";
		cerr << "\t\t[bootstrap opt == 0]: Don't bootstrap\n";
		cerr << "\t\t[bootstrap opt == 1]: Bootstrap\n";
		cerr << "**********************************\n";
		return -1;
	}

	// Random seed for offsets
	srand (time(NULL));
	//srand(1);	

	///////////////////////////////////////////////////////////////////////////
	// Load the data structure that we will use:
	std::ifstream f;
	std::string line;
	f.open("../data/g10_data.txt");
	f.ignore(1000, '\n'); // skip 1 line of the file
	if( f.is_open() ){
		while( getline(f,line) ){
			std::istringstream ss(line);
			string flag;
			double g, alphab, errb, alphat, errt;
			ss >> flag >> g >> alphab >> errb >> alphat >> errt;
			if( flag == "75scut" ){
				DataPoint dat;
				dat.G = g;
				dat.AlphaTop = alphat;
				dat.AlphaBot = alphab;
				dat.AlphaErrTop = errt;
				dat.AlphaErrBot = errb;
				g10_data.push_back( dat );
			}
		}
	}
	f.close();
	nDataPoints = g10_data.size();

	bootstrap_opt = atoi(argv[1]);
	if( bootstrap_opt != 0 && bootstrap_opt != 1 ){ cerr << "unexpected input!\n"; exit(-1); }
	if( bootstrap_opt){
		/////////////////////////////////////////////////////////////////////////////
		// Create RANDOM samples of data sets for boostrapping error estimation
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator

		std::uniform_int_distribution<> distr_g10_data(0, g10_data.size()-1);
		for( int i = 0 ; i < g10_data.size() ; ++i ){ 	
			int rand_elem = distr_g10_data(gen);
			bootstrap_g10_data.push_back( g10_data.at(rand_elem) );
		}
		g10_data.clear();
		g10_data = bootstrap_g10_data;
	}

	///////////////////////////////////////////////////////////////////////////
	// Now we can setup our Chi2 function and minimize it!
	//int nParams = 5; // for energy spectrum & alpha & shift
	//int nParams = 3; // for energy spectrum & alpha & shift
	int nParams = 4; // for energy spectrum & alpha & shift


	double *params = new double[nParams];
	// Set initial guess for params
	//params[0] = 2;		// a
	//params[params[1] = 1.5;	// b
	//params[params[2] = 95;		// Emax
	//params[params[3] = 0.86;	// alpha0
	//params[params[4] = 2;		// pT/cm offset
	//params[0] = 80;		// Emax
	//params[1] = 0.86;	// alpha0
	//params[2] = 2;		// pT/cm offset
	params[0] = 2;		// a
	params[1] = 1.5;	// b
	params[2] = 0.86;	// alpha0
	params[3] = 2;		// pT/cm offset

	// Define some working memory for the minimization algorithm
	nDataPoints+=1;
	double 	*fvec	=	new double[nDataPoints];
	int 	*iwork	=	new int[nDataPoints];
	int 	sdwork	=	nDataPoints*nParams + 5*nParams+nDataPoints;
	double 	*dwork	=	new double[sdwork];

	///////////////////////////////////////////////////////////////////////////
	// Minimize!
	num_pars = nParams;
	

	//cout << "...minimizing...\n";
	int r	=	lmdif1(Chi2,	NULL,
			nDataPoints,	nParams,	params,
			fvec,	1e-7,
			iwork,	dwork,	sdwork);
	//cout << "...done minimizing\n";
	//cout << "r value of minimization: " << r << "\n";
	// Print out the final parameters found:
	cout << "\n";
	cout << "MinPars:\t";
	for( int i = 0 ; i < nParams ; ++i ){
		cout << std::setprecision(15); 
		cout << params[i] << "\t";
	}
	cout << "\n";
	for( double G = -300 ; G <= 300 ; G += 1 ){
		double theo = alpha(G,params);
		cout << G << " " << theo << "\n";
	}

	cout << "\n";



	///////////////////////////////////////////////////////////////////////////
	// Cleanup
	delete []params;
	delete []fvec;
	delete []iwork;
	delete []dwork;

	return 1;
}

int Chi2( void *p, int m, int n, const double *pars, double *fvec, int iflag){
	// p not used??
	// m = num data points
	// n = 0, not used??
	// x = parameters
	// fvec = residuals
	// iflag not used??

	// Want to return chi2 which is calculated given the parameters to be minimized
	std::cerr << "************ Fitting in progress *************\n";
	int dat_point = 0;

	///////////////////////////////////////////////////////////////////////////
	// Perform loop over g10 data
	double chi2 = 0.;
	for( int i = 0 ; i < g10_data.size(); ++i ){
		double dAlpha = g10_data.at(i).AlphaBot;
		double dG = g10_data.at(i).G;
		double dAlphaErr = g10_data.at(i).AlphaErrBot;


		// Integral over energy spectrum
		double theo = alpha(dG,pars);


		double this_chi2 = (dAlpha - theo) / dAlphaErr;
		fvec[dat_point] = this_chi2;
		++dat_point;
		chi2 += pow(this_chi2,2);
	} // end loop over g10_data


	fvec[dat_point] = 0;
	++dat_point;

	std::cerr << "------------Finished calculations!------------\n";
	std::cerr << "\tOptions used:\n";
	std::cerr << "\t\tgrad_order = " << grad_order << "\n";
	std::cerr << "\tCurrent parameters: " << num_pars << "\n";
	for( int i = 0 ; i < num_pars ; ++i ) std::cerr << "\t\t" << pars[i] << "\n";
	std::cerr << "\tCurrent chi2: " << chi2 << "\n";
	std::cerr << "\tReduced chi2: " << chi2 / (nDataPoints-num_pars) << "\n";
	std::cerr << "\tCurrent residuals: \n";
	for(int i = 0; i < m ; ++i) std::cerr << "\t\t" << fvec[i] << "\n";
	std::cerr << "**********************************************\n\n";

	return 0.;
}

double espec( double e, double a, double b, double emax ){
	return pow(e,a)*pow(emax-e,b);
}

double zbar( double en ){
	const double H = 12.;				// cm

	double z = 0;
	// delta-z calculation:
	if( (en / 1.02) < H ){
		z = 2./5. * (en/1.02) - H/2.;
	}
	else{
		double eta = pow(1. - H/ (en/1.02) , 3./2.);
		z = (en/1.02)/(1. - eta)*(2./5. - eta + 3./5.*pow(eta,5./3.) ) - H/2.;
	}

	return z; // in cm
}

double alpha( double G, const double *pars ){
	const double en_step = 0.1; 			// neV
	const double gamman = 29.1647E-6*2*M_PI;	// Hz/pT
	const double T = 180.; 				// s
	const double e_min = 0*1.02;			// 12cm -> neV

	double norm = 0;
	//double a = pars[0];
	//double b = pars[1];
	//double emax = pars[2];
	//double alpha0 = pars[3];
	//double shift = pars[4];
	//double a = 2;
	//double b = 1.5;
	//double emax = pars[0];
	//double alpha0 = pars[1];
	//double shift = pars[2];
	double a = pars[0];
	double b = pars[1];
	double emax = 95;
	double alpha0 = pars[2];
	double shift = pars[3];

	// First we need to calculate <z>:
	double avg_z = 0;
	for(double en=e_min; en <= emax ; en += en_step ){
		// energy spectrum normalization:
		norm += ( espec(en,a,b,emax) * en_step );
		
		avg_z += (zbar(en) * espec(en,a,b,emax) * en_step);
	}
	avg_z /= norm;


	// Now we can calculate alpha(G):
	double theo = 0;
	for(double en=e_min; en <= emax ; en += en_step ){
		double deltaZ = zbar(en) - avg_z;
		theo += ( cos( gamman * (G-shift) * T * deltaZ ) * espec(en,a,b,emax) * en_step );
	}
	theo *= alpha0;
	theo /= norm;
	return theo;
}
