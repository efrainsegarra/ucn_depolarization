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
double alpha_g10( double G, const double *pars );
double alpha_g11( double G, const double *pars );
double alpha_g1m1( double G, const double *pars );
double alpha_g20( double G, const double *pars );
double zbar( double en );
double espec( double e, double a, double b, double max );
double tau_11(double en, double p_electrode, double p_ring);
double tau_1m1(double en, double p_electrode, double p_ring);
double tau_20(double en, double p_electrode, double p_ring);
int Chi2( void *p, int m, int n, const double *pars, double *fvec, int iflag);

int nDataPoints;
int num_pars = -1;
int bootstrap_opt = -1;
int bottop_opt = -1;
int diffuse_opt = -1;
struct DataPoint{
	double G;
	double AlphaTop;
	double AlphaErrTop;
	double AlphaBot;
	double AlphaErrBot;
};
std::vector<DataPoint> g10_data;
std::vector<DataPoint> g11_data;
std::vector<DataPoint> g1m1_data;
std::vector<DataPoint> g20_data;
std::vector<DataPoint> bootstrap_g10_data;
std::vector<DataPoint> bootstrap_g11_data;
std::vector<DataPoint> bootstrap_g1m1_data;
std::vector<DataPoint> bootstrap_g20_data;

void load_data();
const double en_step = 0.1; 			// neV
const double gamman = 29.1647E-6*2*M_PI;	// Hz/pT
const double T = 180.; 				// s
const double R = 40;				// cm
const double H = 18;				// cm
const double mN = 1.67492750056E-27;		// kg
const double emax = 54;				// neV

// Main
int main(int argc, char ** argv){

	if (argc<4){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [bootstrap opt] [bot_top opt] [diffuse opt]\n";
		cerr << "**********************************\n";
		cerr << "\t\t[bootstrap opt == 0]: Don't bootstrap\n";
		cerr << "\t\t[bootstrap opt == 1]: Bootstrap\n";
		cerr << "**********************************\n";
		cerr << "\t\t[bot_top opt == 0]: Fit bottom chamber\n";
		cerr << "\t\t[bot_top opt == 1]: Fit top chamber\n";
		cerr << "**********************************\n";
		cerr << "\t\t[diffuse opt == 0]: p_ring=p_electrode\n";
		cerr << "\t\t[diffuse opt == 1]: p_ring=0\n";
		cerr << "\t\t[diffuse opt == 2]: p_ring=1\n";
		cerr << "**********************************\n";
		return -1;
	}
	bootstrap_opt = atoi(argv[1]);
	bottop_opt = atoi(argv[2]);
	diffuse_opt = atoi(argv[3]);
	if( bootstrap_opt != 0 && bootstrap_opt != 1 )		{ cerr << "unexpected input!\n"; exit(-1); }
	if( bottop_opt != 0 && bottop_opt != 1 )		{ cerr << "unexpected input!\n"; exit(-1); }
	if( diffuse_opt < 0 || diffuse_opt > 2 )		{ cerr << "unexpected input!\n"; exit(-1); }

	// Random seed for offsets
	srand (time(NULL));
	//srand(1);	

	///////////////////////////////////////////////////////////////////////////
	// Load the data structure that we will use:
	load_data();
	nDataPoints = g10_data.size() + g1m1_data.size() + g11_data.size() + g20_data.size();


	///////////////////////////////////////////////////////////////////////////
	// Now we can setup our Chi2 function and minimize it!
	//int nParams = 5; // for energy spectrum & alpha & shift
	//int nParams = 3; // for energy spectrum & alpha & shift
	int nParams = 8; // for energy spectrum & alpha & shifts & diffusivities


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
	params[2] = 1;		// alpha0
	params[3] = 0;		// g10 pT/cm offset
	params[4] = 0;		// g1m1 pT/cm offset
	params[5] = 0;		// g11 pT/cm offset
	params[6] = 0;		// g20 pT/cm offset
	params[7] = 0.5;	// diffusivity
	//params[8] = 0.5;	// diffusivity

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
	//for( double G = -300 ; G <= 300 ; G += 1 ){
	//for( do	double theo = alpha(G,params);
	//for( do	cout << G << " " << theo << "\n";
	//for( do}

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
	double chi2_g10 = 0.;
	for( int i = 0 ; i < g10_data.size(); ++i ){
		double dAlpha, dAlphaErr;
		double dG 		= g10_data.at(i).G;
		if( bottop_opt == 0 ){
			dAlpha 		= g10_data.at(i).AlphaBot;
			dAlphaErr 	= g10_data.at(i).AlphaErrBot;
		}
		else{
			dAlpha 		= g10_data.at(i).AlphaTop;
			dAlphaErr 	= g10_data.at(i).AlphaErrTop;
		}

		// Integral over energy spectrum
		double theo = alpha_g10(dG,pars);


		double this_chi2 = (dAlpha - theo) / dAlphaErr;
		cout << "\tG10 " << dG << " " << dAlpha << " " << theo << " " << this_chi2 << "\n";
		fvec[dat_point] = this_chi2;
		++dat_point;
		chi2_g10 += pow(this_chi2,2);
	} // end loop over g10_data

	///////////////////////////////////////////////////////////////////////////
	// Perform loop over g1m1 data
	double chi2_g1m1 = 0.;
	for( int i = 0 ; i < g1m1_data.size(); ++i ){
		double dAlpha, dAlphaErr;
		double dG 		= g1m1_data.at(i).G;
		if( bottop_opt == 0 ){
			dAlpha 		= g1m1_data.at(i).AlphaBot;
			dAlphaErr 	= g1m1_data.at(i).AlphaErrBot;
		}
		else{
			dAlpha 		= g1m1_data.at(i).AlphaTop;
			dAlphaErr 	= g1m1_data.at(i).AlphaErrTop;
		}

		// Integral over energy spectrum
		double theo = alpha_g1m1(dG,pars);


		double this_chi2 = (dAlpha - theo) / dAlphaErr;
		cout << "\tG1m1 " << dG << " " << dAlpha << " " << theo << " " << this_chi2 << "\n";
		fvec[dat_point] = this_chi2;
		++dat_point;
		chi2_g1m1 += pow(this_chi2,2);
	} // end loop over g1m1_data

	///////////////////////////////////////////////////////////////////////////
	// Perform loop over g11 data
	double chi2_g11 = 0.;
	for( int i = 0 ; i < g11_data.size(); ++i ){
		double dAlpha, dAlphaErr;
		double dG 		= g11_data.at(i).G;
		if( bottop_opt == 0 ){
			dAlpha 		= g11_data.at(i).AlphaBot;
			dAlphaErr 	= g11_data.at(i).AlphaErrBot;
		}
		else{
			dAlpha 		= g11_data.at(i).AlphaTop;
			dAlphaErr 	= g11_data.at(i).AlphaErrTop;
		}

		// Integral over energy spectrum
		double theo = alpha_g11(dG,pars);


		double this_chi2 = (dAlpha - theo) / dAlphaErr;
		cout << "\tG11 " << dG << " " << dAlpha << " " << theo << " " << this_chi2 << "\n";
		fvec[dat_point] = this_chi2;
		++dat_point;
		chi2_g11 += pow(this_chi2,2);
	} // end loop over g11_data


	///////////////////////////////////////////////////////////////////////////
	// Perform loop over g20 data
	double chi2_g20 = 0.;
	for( int i = 0 ; i < g20_data.size(); ++i ){
		double dAlpha, dAlphaErr;
		double dG 		= g20_data.at(i).G;
		if( bottop_opt == 0 ){
			dAlpha 		= g20_data.at(i).AlphaBot;
			dAlphaErr 	= g20_data.at(i).AlphaErrBot;
		}
		else{
			dAlpha 		= g20_data.at(i).AlphaTop;
			dAlphaErr 	= g20_data.at(i).AlphaErrTop;
		}

		// Integral over energy spectrum
		double theo = alpha_g20(dG,pars);


		double this_chi2 = (dAlpha - theo) / dAlphaErr;
		cout << "\tG20 " << dG << " " << dAlpha << " " << theo << " " << this_chi2 << "\n";
		fvec[dat_point] = this_chi2;
		++dat_point;
		chi2_g20 += pow(this_chi2,2);
	} // end loop over g20_data


	fvec[dat_point] = 0;
	++dat_point;
	
	std::cerr << "------------Finished calculations!------------\n";
	std::cerr << "\tOptions used:\n";
	std::cerr << "\t\tbootstrap_opt = " << bootstrap_opt << "\n";
	std::cerr << "\t\tbottop_opt = " << bottop_opt << "\n";
	std::cerr << "\tCurrent parameters: " << num_pars << "\n";
	for( int i = 0 ; i < num_pars ; ++i ) std::cerr << "\t\t" << pars[i] << "\n";
	double chi2 = chi2_g10 + chi2_g11 + chi2_g1m1 + chi2_g20;
	std::cerr << "\tTotal chi2: " << chi2 << "\n";
	std::cerr << "\tReduced chi2: " << chi2 / (nDataPoints-num_pars) << "\n";
	std::cerr << "\tIndividual chi2:\n";
	std::cerr << "\t\tG10: " << chi2_g10 << "\t" << chi2_g10/(g10_data.size()-4) << "\n";
	std::cerr << "\t\tG1m1: " << chi2_g1m1 << "\t" << chi2_g1m1/(g1m1_data.size()-5) << "\n";
	std::cerr << "\t\tG11: " << chi2_g11 << "\t" << chi2_g11/(g11_data.size()-5) << "\n";
	std::cerr << "\t\tG20: " << chi2_g20 << "\t" << chi2_g20/(g20_data.size()-5) << "\n";
	std::cerr << "\tCurrent residuals: \n";
	for(int i = 0; i < m ; ++i) std::cerr << "\t\t" << fvec[i] << "\n";
	std::cerr << "**********************************************\n\n";

	return 0.;
}

double espec( double e, double a, double b, double max ){
	return pow(e,a)*pow(max-e,b);
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

double alpha_g10( double G, const double *pars ){
	double e_min;
	if( bottop_opt == 0 ){
		e_min = 0*1.02;
	}
	else{
		e_min = 12*1.02;
	}

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
	double alpha0 = pars[2];
	double g10_shift = pars[3];
	double g1m1_shift = pars[4];
	double g11_shift = pars[5];
	double g20_shift = pars[6];
	//double p_electrode = pars[7];
	//double p_ring = pars[8];
	double p = pars[7];
	double p_electrode, p_ring;
	if( diffuse_opt == 0 ){
		p_electrode = p;
		p_ring = p;
	}
	else if(diffuse_opt == 1 ){
		p_electrode = p;
		p_ring = 0;
	}
	else if(diffuse_opt == 2 ){
		p_electrode = p;
		p_ring = 1;
	}

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
		theo += ( cos( gamman * (G-g10_shift) * T * deltaZ ) * espec(en,a,b,emax) * en_step );
	}
	theo *= alpha0;
	theo /= norm;
	return theo;
}

double tau_11(double en, double p_electrode, double p_ring){
	double L = 10;	// cm
	double e1 = 63;	// cm
	double r1 = 24;	// cm
	// convert energy (neV) to J
	double v = sqrt( 2 * en*1.60217733000001E-28 / mN ) * 100; // m/s --> cm/s
	return ( L + sqrt( pow(e1*p_electrode,2) + pow(r1*p_ring,2) ) )/v;
}
double tau_1m1(double en, double p_electrode, double p_ring){
	return tau_11(en,p_electrode,p_ring);
}
double tau_20(double en, double p_electrode, double p_ring){
	double e2 = 22;		// cm
	double e2p = 15;	// cm
	double r2p = 9;		// cm
	// convert energy (neV) to J
	double v = sqrt( 2 * en*1.60217733000001E-28 / mN ) * 100; // m/s --> cm/s
	return ( e2*p_electrode + pow(pow(p_electrode/e2p,2)+pow(p_ring/r2p,2),-0.5))/v;
}	

double alpha_g11( double G, const double *pars ){
	double e_min;
	if( bottop_opt == 0 ){
		e_min = 0*1.02;
	}
	else{
		e_min = 12*1.02;
	}

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
	double alpha0 = pars[2];
	double g10_shift = pars[3];
	double g1m1_shift = pars[4];
	double g11_shift = pars[5];
	double g20_shift = pars[6];
	//double p_electrode = pars[7];
	//double p_ring = pars[8];
	double p = pars[7];
	double p_electrode, p_ring;
	if( diffuse_opt == 0 ){
		p_electrode = p;
		p_ring = p;
	}
	else if(diffuse_opt == 1 ){
		p_electrode = p;
		p_ring = 0;
	}
	else if(diffuse_opt == 2 ){
		p_electrode = p;
		p_ring = 1;
	}

	// First we need to calculate <z>:
	for(double en=e_min; en <= emax ; en += en_step ){
		// energy spectrum normalization:
		norm += ( espec(en,a,b,emax) * en_step );
	}


	// Now we can calculate alpha(G):
	double theo = 0;
	for(double en=e_min; en <= emax ; en += en_step ){
		if( G == 0 ){
			theo += ( espec(en,a,b,emax) * en_step );
		}
		else{
			theo += ( exp( - gamman*gamman * R*R * T * pow(G-g11_shift,2) * tau_11(en,p_electrode,p_ring) / 4. ) * espec(en,a,b,emax) * en_step );
		}
	}
	theo *= alpha0;
	theo /= norm;
	return theo;
}
double alpha_g1m1( double G, const double *pars ){
	double e_min;
	if( bottop_opt == 0 ){
		e_min = 0*1.02;
	}
	else{
		e_min = 12*1.02;
	}

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
	double alpha0 = pars[2];
	double g10_shift = pars[3];
	double g1m1_shift = pars[4];
	double g11_shift = pars[5];
	double g20_shift = pars[6];
	//double p_electrode = pars[7];
	//double p_ring = pars[8];
	double p = pars[7];
	double p_electrode, p_ring;
	if( diffuse_opt == 0 ){
		p_electrode = p;
		p_ring = p;
	}
	else if(diffuse_opt == 1 ){
		p_electrode = p;
		p_ring = 0;
	}
	else if(diffuse_opt == 2 ){
		p_electrode = p;
		p_ring = 1;
	}

	// First we need to calculate <z>:
	for(double en=e_min; en <= emax ; en += en_step ){
		// energy spectrum normalization:
		norm += ( espec(en,a,b,emax) * en_step );
	}


	// Now we can calculate alpha(G):
	double theo = 0;
	for(double en=e_min; en <= emax ; en += en_step ){
		if( G == 0 ){
			theo += ( espec(en,a,b,emax) * en_step );
		}
		else{
			theo += ( exp( - gamman*gamman * R*R * T * pow(G-g1m1_shift,2) * tau_1m1(en,p_electrode,p_ring) / 4. ) * espec(en,a,b,emax) * en_step );
		}
	}
	theo *= alpha0;
	theo /= norm;
	return theo;
}


double alpha_g20( double G, const double *pars ){
	double e_min;
	if( bottop_opt == 0 ){
		e_min = 0*1.02;
	}
	else{
		e_min = 12*1.02;
	}

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
	double alpha0 = pars[2];
	double g10_shift = pars[3];
	double g1m1_shift = pars[4];
	double g11_shift = pars[5];
	double g20_shift = pars[6];
	//double p_electrode = pars[7];
	//double p_ring = pars[8];
	double p = pars[7];
	double p_electrode, p_ring;
	if( diffuse_opt == 0 ){
		p_electrode = p;
		p_ring = p;
	}
	else if(diffuse_opt == 1 ){
		p_electrode = p;
		p_ring = 0;
	}
	else if(diffuse_opt == 2 ){
		p_electrode = p;
		p_ring = 1;
	}

	// First we need to calculate <z>:
	for(double en=e_min; en <= emax ; en += en_step ){
		// energy spectrum normalization:
		norm += ( espec(en,a,b,emax) * en_step );
	}


	// Now we can calculate alpha(G):
	double theo = 0;
	for(double en=e_min; en <= emax ; en += en_step ){
		if( G == 0 ){
			theo += ( espec(en,a,b,emax) * en_step );
		}
		else{
			theo += ( exp( - gamman*gamman * (4./45.*pow(H,4) + 1./48.*pow(R,4) ) * T * pow(G-g20_shift,2) * tau_20(en,p_electrode,p_ring) ) * espec(en,a,b,emax) * en_step );
		}
	}
	theo *= alpha0;
	theo /= norm;
	return theo;
}


void load_data(){
	std::ifstream f;
	std::string line;

	// G10
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
				dat.AlphaErrTop = errt*1;
				dat.AlphaErrBot = errb*1;
				g10_data.push_back( dat );
			}
		}
	}
	f.close();

	// G1m1
	f.open("../data/g1m1_data.txt");
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
				dat.AlphaErrTop = errt*1;
				dat.AlphaErrBot = errb*1;
				g1m1_data.push_back( dat );
			}
		}
	}
	f.close();

	// G11
	f.open("../data/g11_data.txt");
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
				dat.AlphaErrTop = errt*1;
				dat.AlphaErrBot = errb*1;
				g11_data.push_back( dat );
			}
		}
	}
	f.close();

	// G20
	f.open("../data/g20_data.txt");
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
				dat.AlphaErrTop = errt*1;
				dat.AlphaErrBot = errb*1;
				g20_data.push_back( dat );
			}
		}
	}
	f.close();


	
	/////////////////////////////////////////////////////////////////////////////
	// Create RANDOM samples of data sets for boostrapping error estimation
	if( bootstrap_opt){
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator

		//////////////////////////////////////////////////////////////////////

		std::uniform_int_distribution<> distr_g10_data(0, g10_data.size()-1);
		for( int i = 0 ; i < g10_data.size() ; ++i ){ 	
			int rand_elem = distr_g10_data(gen);
			bootstrap_g10_data.push_back( g10_data.at(rand_elem) );
		}
		g10_data.clear();
		g10_data = bootstrap_g10_data;

		//////////////////////////////////////////////////////////////////////

		std::uniform_int_distribution<> distr_g1m1_data(0, g1m1_data.size()-1);
		for( int i = 0 ; i < g1m1_data.size() ; ++i ){ 	
			int rand_elem = distr_g1m1_data(gen);
			bootstrap_g1m1_data.push_back( g1m1_data.at(rand_elem) );
		}
		g1m1_data.clear();
		g1m1_data = bootstrap_g1m1_data;

		//////////////////////////////////////////////////////////////////////

		std::uniform_int_distribution<> distr_g11_data(0, g11_data.size()-1);
		for( int i = 0 ; i < g11_data.size() ; ++i ){ 	
			int rand_elem = distr_g11_data(gen);
			bootstrap_g11_data.push_back( g11_data.at(rand_elem) );
		}
		g11_data.clear();
		g11_data = bootstrap_g11_data;

		//////////////////////////////////////////////////////////////////////

		std::uniform_int_distribution<> distr_g20_data(0, g20_data.size()-1);
		for( int i = 0 ; i < g20_data.size() ; ++i ){ 	
			int rand_elem = distr_g20_data(gen);
			bootstrap_g20_data.push_back( g20_data.at(rand_elem) );
		}
		g20_data.clear();
		g20_data = bootstrap_g20_data;

	}
}
