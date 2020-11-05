// ========================================================================================
// Perform a Bayesian fit to a set of data points with a generic model
// ========================================================================================

#include "Cosmology.h"
#include "Data1D.h"
#include "Posterior.h"
#include "../FitDeltaSigma/nfwLens.h"
#include "TaperedCovarianceMatrix.h"

using namespace std;

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;

// ... my pars
// std:: string method_Pk;// = "EisensteinHu"; // testing HU, CAMB is slower ... this is only linear
// std:: string method_Pk = "CAMB";
// bool store = false; 
// bool do_nonlinear = true;
// std:: string tt;
const double tiny =1e-4;

int main () {
  try {
    std:: cout << " Save the data model from chain " << std:: endl;
    
    /***** COSMOLOGY *****/  /* BEGIN */
    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck18_};    
    // as used to comupute the stacked shear profiles
    cosmology.set_Omega(0.3);
    cosmology.set_H0(70.0);
    /***** COSMOLOGY *****/  /* END */
    
    // print cosmological parameters
    cosmology.print_parameters();
    
    std:: string filin;
    std:: cin >> filin;
    double zl;
    std:: cin >> zl;
    std:: ifstream ifilin;
    
    double med_lm200, q18_lm200, q82_lm200;
    double med_c200, q18_c200, q82_c200;
    double med_f_off, q18_f_off, q82_f_off;
    double med_sigma_off, q18_sigma_off, q82_sigma_off;
    double med_omega_m, q18_omega_m, q82_omega_m;
    
    int numLines = 0;
    ifstream in(filin.c_str());
    std::string unused;
    while ( std::getline(in, unused) )
      ++numLines;
    std:: cout <<  " number of lines in the chain parameters file " << numLines << std :: endl;
    
    if(numLines==5){
      ifilin.open(filin.c_str());
      if(ifilin.is_open()){
	std:: string line;      
	int nlines_to_skip = 1;
	int nl=0;
	for(int i=0;i<nlines_to_skip;i++){
	  getline(ifilin,line);
	  nl++;
	}
	std:: string sbut;
	double dbut;
	
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_lm200 >> q18_lm200 >> q82_lm200;
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_c200 >> q18_c200 >> q82_c200;
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_f_off >> q18_f_off >> q82_f_off;
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_sigma_off >> q18_sigma_off >> q82_sigma_off;
	
	std:: cout << " lm200     = " << med_lm200 << "  "  << q18_lm200 << "  " << q82_lm200 << std:: endl;
	std:: cout << " c200      = " << med_c200 << "  " << q18_c200 << "  " << q82_c200 << std:: endl;
	std:: cout << " f_off     = " << med_f_off << "  " << q18_f_off << "  " << q82_f_off << std:: endl;
	std:: cout << " sigma_off = " << med_sigma_off << "  " << q18_sigma_off << "  " << q82_sigma_off << std:: endl;
	
      }else{
	std:: cout << filin << " does not exist " << std:: endl;
	std:: cout << " I will STOP here!!!" << std:: endl;
	exit(1);
      }
    }
    if(numLines==6){
      ifilin.open(filin.c_str());
      if(ifilin.is_open()){
	std:: string line;      
	int nlines_to_skip = 1;
	int nl=0;
	for(int i=0;i<nlines_to_skip;i++){
	  getline(ifilin,line);
	  nl++;
	}
	std:: string sbut;
	double dbut;
	
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_lm200 >> q18_lm200 >> q82_lm200;
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_c200 >> q18_c200 >> q82_c200;
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_f_off >> q18_f_off >> q82_f_off;
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_sigma_off >> q18_sigma_off >> q82_sigma_off;
	ifilin >> sbut >> sbut >> dbut >> dbut >> med_omega_m >> q18_omega_m >> q82_omega_m;                              	
	
	std:: cout << " lm200     = " << med_lm200 << "  "  << q18_lm200 << "  " << q82_lm200 << std:: endl;
	std:: cout << " c200      = " << med_c200 << "  " << q18_c200 << "  " << q82_c200 << std:: endl;
	std:: cout << " f_off     = " << med_f_off << "  " << q18_f_off << "  " << q82_f_off << std:: endl;
	std:: cout << " sigma_off = " << med_sigma_off << "  " << q18_sigma_off << "  " << q82_sigma_off << std:: endl;
	std:: cout << "  omega_m  = " << med_omega_m << "  " << q18_omega_m << "  " << q82_omega_m << std:: endl;
	
	cosmology.set_Omega(med_omega_m);
	// print cosmological parameters
	cosmology.print_parameters();
	
      }else{
	std:: cout << filin << " does not exist " << std:: endl;
	std:: cout << " I will STOP here!!!" << std:: endl;
	exit(1);
      }
    }
    
    // halo mass, redshift and concentration
    double m200 = pow(10.0,med_lm200);
    double redshift = zl;
    double c200 = med_c200;
    
    // compute the critical density at redshift   
    double H = cosmology.HH(redshift)/cosmology.HH(0.)*100./3.0857e+19; // in sec^-1
    double rho_crit = 3*H*H/8/M_PI/6.6732e-8/1.98892e+33*3.0857e+24*3.0857e+24*3.0857e+24; // in h^2*M_sun/Mpc^3
    std:: cout << "critical density (z=" << redshift << ") = " <<  rho_crit << std:: endl;
    std:: cout << "  " << std:: endl;
    
    // compute the 3D density profile
    double logr_min = -1.;
    double logr_max = 1.5;
    int step = 128;
    
    std:: cout << redshift << "  " << m200 << "  " << c200 << std:: endl;
    nfwLens lens (&cosmology, redshift, m200, c200, 3.0, med_f_off, med_sigma_off, true);
    if(numLines==6) cosmology.set_Omega(med_omega_m - q18_omega_m);    
    nfwLens lens0(&cosmology, redshift, pow(10.,med_lm200-q18_lm200), med_c200-q18_c200, 3.0, med_f_off-q18_f_off, med_sigma_off-q18_sigma_off, true);
    if(numLines==6) cosmology.set_Omega(med_omega_m + q82_omega_m);        
    nfwLens lens1(&cosmology, redshift, pow(10.,med_lm200+q82_lm200), med_c200+q82_c200, 3.0, med_f_off+q82_f_off, med_sigma_off+q82_sigma_off, true);
    
    std:: cout << " running " << std:: endl;
    
    std::vector<double> lr = cbl::linear_bin_vector(step,logr_min,logr_max);	
    std::vector<double> r(step);    
    for(int i=0;i<step;i++){
      r[i] = pow(10.,lr[i]);
      std:: cout << r[i] << "   " 
		 << lens.deltasigma_1h(r[i]) << "  " << lens.deltasigma_2h(r[i]) << "  "
		 << lens.deltasigma(r[i]) << "  "
	
		 << lens0.deltasigma_1h(r[i]) << "  " << lens0.deltasigma_2h(r[i]) << "  "
		 << lens0.deltasigma(r[i]) << "  "
	
		 << lens1.deltasigma_1h(r[i]) << "  " << lens1.deltasigma_2h(r[i]) << "  "
		 << lens1.deltasigma(r[i]) << "  "
	
		 << std:: endl;
    }
  }
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  std:: cout << " ... end of work ... ;-) " << std:: endl;
  return 0;  
}  






