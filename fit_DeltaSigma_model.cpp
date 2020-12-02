// ========================================================================================
// Perform a Bayesian fit to a set of data points with a generic model
// ========================================================================================

#include "Cosmology.h"
#include "Data1D.h"
#include "Posterior.h"
#include "../FitDeltaSigma/nfwLens.h"
#include "TaperedCovarianceMatrix.h"

template <class T>
int locate (const std::vector<T> &v, const T x){
  size_t n = v.size ();
  int jl = -1;
  int ju = n;
  bool as = (v[n-1] >= v[0]);
  while (ju-jl > 1){
    int jm = (ju+jl)/2;
    if ((x >= v[jm]) == as)
      jl=jm;
    else
      ju=jm;
  }
  if (x == v[0])
    return 0;
  else if (x == v[n-1])
    return n-2;
  else
    return jl;
}

double getY(std:: vector<double> x, std:: vector<double> y,double xi){
  int nn = x.size();                                                               
  if(x[0]<x[nn-1]){                                                                                     
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }      
  else{                                                                                                            
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];        
  }                                    
  int i = locate (x,xi);           
  i = std::min (std::max (i,0), int (nn)-2);
  double f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    double a0,a1,a2,a3,f2;                                                       
    f2 = f*f;                                                    
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];            
    a1 = y[i-1] - y[i] - a0;                                                                              
    a2 = y[i+1] - y[i-1];                                                              
    a3 = y[i];                                                                                     
    return a0*f*f2+a1*f2+a2*f+a3;                                                    
  }                                                                      
  else return f*y[i+1]+(1-f)*y[i];           
}

using namespace std;

std:: vector<double> lrr, D2h;

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

const double min_ngal_in_bin = 10;
const double rmin = 0.01;

int main () {
  try {
    std:: cout << " Save the data model from chain " << std:: endl;
    
    /***** COSMOLOGY *****/  /* BEGIN */
    // cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck18_};    
    // as used to comupute the stacked shear profiles
    // cosmology.set_Omega(0.3);
    // cosmology.set_H0(70.0);
    /*************************************************************************************************************************/
    /*************************************************************************************************************************/
    // 300th Cosmology
    double Omega_matter = 0.307;
    double Omega_baryon = 0.048;
    const double Omega_nu = 0.;
    double massless_neutrinos = 2.04;
    const double Omega_radiation = 0.;
    int massive_neutrinos = 1.0;
    double Omega_DE = 0.693;
    double hh = 0.678;
    double scalar_amp = 2.139e-9;
    double scalar_pivot = 0.05;
    double n_spec = 0.96;
    // double tau = 0.066;
    const double wa = 0.;
    const double w0 = -1.;
    //
    /***** COSMOLOGY *****/  /* BEGIN */
    // cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};
    cbl::cosmology::Cosmology cosmology {Omega_matter, Omega_baryon,
                                        Omega_nu, massless_neutrinos, massive_neutrinos, Omega_DE, Omega_radiation,
                                        hh, scalar_amp, scalar_pivot, n_spec, w0, wa};
    // as used to comupute the stacked shear profiles
    std:: cout << "s8 EH " << cosmology.sigma8_Pk("EisensteinHu",0.0) << std:: endl;       
    /***** COSMOLOGY *****/  /* END */
    
    // print cosmological parameters
    cosmology.print_parameters();
    
    std:: string filin;
    std:: cin >> filin;
    double zl;
    std:: cin >> zl;
    
    // soon after it checks if there is a tabulated file for the 2halo term
    std:: ostringstream osz;
    osz << zl;    
    std:: string sz;
    sz = osz.str();
    std:: string filin_table = "table_2h_model_z=" + sz + ".txt";
    std:: ifstream ifilin_table;
    ifilin_table.open(filin_table.c_str());
    if(ifilin_table.is_open()){
      std:: cout << " table file for the 2halo term  exists: " << filin_table << std:: endl;
      std:: cout << " I read it!!! " << std:: endl;
      double ri, Di2h;
      while(ifilin_table >> ri >> Di2h){
	lrr.push_back(log10(ri));
	D2h.push_back(Di2h);	      
      }
      std:: cout << " the file containts " << lrr.size() << " lines " << std:: endl;
    }		
    
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

    if(numLines==3){
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
	
	std:: cout << " lm200     = " << med_lm200 << "  "  << q18_lm200 << "  " << q82_lm200 << std:: endl;
	std:: cout << " c200      = " << med_c200 << "  " << q18_c200 << "  " << q82_c200 << std:: endl;

	med_f_off=q18_f_off=q82_f_off=0;
	med_sigma_off=q18_sigma_off=q82_sigma_off=0;
	
      }else{
	std:: cout << filin << " does not exist " << std:: endl;
	std:: cout << " I will STOP here!!!" << std:: endl;
	exit(1);
      }
    }       
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


    std:: string filout = filin + "_median_model.txt";
    std:: ofstream ofilout(filout.c_str());            
    
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
    
    if (numLines>3){
      
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
    }else{
      std:: vector<double> radius;      
      std:: string fileprof;
      cin >> fileprof;
      ifstream ifileprof(fileprof.c_str());
      if(ifileprof.is_open()){
	while (true){
	  double a3, but, a4, a5, ngal_in_bin;
	  ifileprof >> a3  >> but >> but >> but >> but >> but >> a4 >> a5 >> ngal_in_bin;	 
	  if(ifileprof.eof()) break;
	  if (a3 > rmin){
	    // in this case we do not account for the SL region in simulations
	    //if(a4 > 0){
	    if(ngal_in_bin >= min_ngal_in_bin){
	      radius.push_back(a3);
	      // std::cout << a3 << "  " << a4 << "  " << a5 << std::endl;
	      //}
	    }
	  }
	}
	ifileprof.close();	
      }else{
	std:: cout << fileprof << " does not exist ... check this out!!! " << std:: endl;
	exit(1);
      }
      int step = radius.size();
      if(lrr.size()>0){
	double halo_bias = cosmology.bias_halo(m200,redshift,"Tinker","EisensteinHu");    	
	double halo_bias0 = cosmology.bias_halo(pow(10.,med_lm200-q18_lm200),redshift,"Tinker","EisensteinHu");
	double halo_bias1 = cosmology.bias_halo(pow(10.,med_lm200+q82_lm200),redshift,"Tinker","EisensteinHu");    	
	for(int i=0;i<step;i++){
	  double D2hi = getY(lrr,D2h,log10(radius[i]));	  
	  ofilout << radius[i] << "   "
		  << lens.deltasigma_1h(radius[i])+halo_bias*D2hi << "  " 	  
		  << lens0.deltasigma_1h(radius[i])+halo_bias0*D2hi << "  " 
		  << lens1.deltasigma_1h(radius[i])+halo_bias1*D2hi
		  << std:: endl;
	}
      }else{
	for(int i=0;i<step;i++){
	  ofilout << radius[i] << "   " 
		  << lens.deltasigma_1h(radius[i]) << "  " 	  
		  << lens0.deltasigma_1h(radius[i]) << "  " 
		  << lens1.deltasigma_1h(radius[i]) 
		  << std:: endl;
	}
      }
      ofilout.close();
    }
  }
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  std:: cout << " ... end of work ... ;-) " << std:: endl;
  return 0;  
}  






