//KLASS
#include "ClassEngine.hh"

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

using namespace std;

// example run: one can specify as a second argument a preicison file
int main(int argc, char** argv) {
  
  const int l_max_scalars =  2508;

  ClassParams pars;
  
  pars.add("omega_b", 0.022252);
  pars.add("omega_cdm", 0.11987);
  pars.add("100*theta_s", 1.040778);
  pars.add("tau_reio", 0.0789);
  pars.add("ln10^{10}A_s", 3.0929);
  pars.add("n_s", 0.96475);

  pars.add("annihilation", 0.);

  pars.add("N_ur", 2.0328);
  pars.add("N_ncdm", 1);
  pars.add("m_ncdm", 0.06);

  pars.add("input_verbose", 1);
  pars.add("background_verbose", 1);
  pars.add("thermodynamics_verbose", 1);
  pars.add("perturbations_verbose", 1);
  pars.add("transfer_verbose", 1);
  pars.add("primordial_verbose", 1);
  pars.add("spectra_verbose", 1);
  pars.add("nonlinear_verbose", 1);
  pars.add("lensing_verbose", 1);
  pars.add("output_verbose", 1);

  pars.add("output", "tCl,pCl,lCl"); // pol +clphi
  pars.add("lensing", true); // note boolean
  pars.add("l_max_scalars", l_max_scalars);
  pars.add("format", "camb");


  ClassEngine* class_engine(0);

  try {
    class_engine = new ClassEngine(pars);
    cout << "[test_class_cpp] CLASS engine initialisation succeeded" << endl;
  }
  catch (exception const &e) {
    cerr << "[test_class_cpp] CLASS engine initialisation failed:\n" << e.what() << endl;
  }

  cout << "[test_class_cpp] Updating CLASS engine..." << endl;

  std::vector<double> new_pars;

  new_pars.push_back(0.022252);
  new_pars.push_back(0.11987);
  new_pars.push_back(1.040778);
  new_pars.push_back(0.0789);
  new_pars.push_back(3.0929);
  new_pars.push_back(0.96475);
  new_pars.push_back(1e-3);

  if (class_engine->updateParValues(new_pars)) {
    cout << "[test_class_cpp] CLASS engine update succeeded" << endl;
  }
  else {
    cerr << "[test_class_cpp] CLASS engine update failed" << endl;
  }

  cout << "[test_class_cpp] Deleting CLASS engine" << endl;
  delete class_engine;
}
