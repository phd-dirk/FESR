#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <stdexcept>

using std::vector;
using std::cout;
using std::endl;
using std::complex;
using std::string;
using std::runtime_error;

#include "./src/constants.hpp"
// Theoretical Moments
#include "./src/adler_function.hpp"
#include "./src/alpha_s.hpp"
#include "./src/numerics.hpp"

// Chisquared
#include "./src/chisquared.hpp"
#include "./src/s0_sets.hpp"
#include "./src/weights.hpp"

// MINUIT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// Mutlidimensional Root-finding
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

// Config
using json = nlohmann::json;

using ROOT::Math::Minimizer;
using ROOT::Math::Factory;
using ROOT::Math::Functor;

// UBLAS
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
using ublas::matrix;
using ublas::prod;

using std::ifstream;
using std::ofstream;

vector<double> readExpMomentsFromFile() {
  vector<double> expMoms;
  std::ifstream expMomFile;
  expMomFile.open("./data/expMom.dat");
  if(!expMomFile) {
    std::cerr << "Unable to open file expMom.dat";
    exit(1);
  }

  double x;
  while (expMomFile >> x) {
    expMoms.push_back(x);
  }
  expMomFile.close();
  return expMoms;
}


// following https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
bool almostEqualRelative(const double &a, const double &b) {
  // Machine epsilon https://en.wikipedia.org/wiki/Machine_epsilon
  double maxRelDiff = 1.0e-14;//2.22e-16;
  double diff = abs(a - b);
  double largest = (b > a) ? b : a;
  if ( diff <= largest * maxRelDiff )
    return true;
  return false;
}

bool compareVectors(const vector<double> &vec1, const vector<double> &vec2) {
  for(std::vector<double>::size_type i = 0; i != vec1.size(); i++) {
    if(!almostEqualRelative(vec1[i], vec2[i])) {
      cout << vec1[i] << " = " << vec2[i] << endl;
      return false;
    }
  }
  return true;
}

bool compareSquareMatrix(const matrix<double> &m1, const matrix<double> &m2) {
  for (auto i = 0; i < m1.size1(); i++) {
    for (auto j = 0; j < m1.size2(); j++) {
      if (!almostEqualRelative(m1(i, j), m2(i, j))) {
        cout << "at m(" << i << ", " << j << ") \t" << m1(i, j) << " != " << m2(i, j) << endl;
        return false;
      }
    }
  }
  return true;
}

complex<double> testFunction(complex<double> z) {
  return 1.0/z;
}

int dof(const json config) {
  int numS0s = config["parameters"]["s0Set"].size();
  int numVar = 0;
  if (!config["variables"]["astau"]["fixed"]) {
    numVar++;
  }
  if (!config["variables"]["aGGInv"]["fixed"]) {
    numVar++;
  }
  if (!config["variables"]["rhoVpA"]["fixed"]) {
    numVar++;
  }
  if (!config["variables"]["c8VpA"]["fixed"]) {
    numVar++;
  }
  return numS0s-numVar;
}

void writeOutput(const string filePath, const double *variables, const double *errors, const double &chi2, const double &edm, const json config) {
  ofstream file;
  file.open(filePath, std::ios::app);
  file << std::setprecision(15);
  file << config["parameters"]["s0Set"].size();

  // add variables and errors
  for(int i = 0; i < 4; i++) {
    file << "," << variables[i] << "," << errors[i];
  }

  // add chi2 and emd
  file << "," << chi2 << "," << chi2/dof(config) << "," << edm;

  // add used s0s
  vector<double> s0s = config["parameters"]["s0Set"];
  std::stringstream ss;
  ss << std::setprecision(15);
  for(size_t i = 0; i < s0s.size(); ++i) {
    if(i != 0)
      ss << " ";
    ss << s0s[i];
  }
  std::string s0sStr = ss.str();
  cout << s0sStr << endl;
  file << ",[" << ss.str() << "]";

  cout << "weight: " << config["parameters"]["weight"] << endl;
  // add used weight
  file << "," << config["parameters"]["weight"];

  file << endl;
  file.close();
}


// test multidimensional root-finding
struct rparams {
  double a;
  double b;
};
int alpha_f(const gsl_vector *x, void *params, gsl_vector *f) {
  double a = ((rparams *) params)->a;
  double b = ((rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);

  const double y0 = 0.6234358730725664 - (0.2222222222222222*x0)/(pow(x0,2) + pow(x1,2)) -
    0.19753086419753085*log(pow(x0,2) + pow(x1,2)) + 0.19753086419753085*log(pow(9. + 16.*x0,2) + 256.*pow(x1,2));
  const double y1 = 0.7853981633974483 + (0.2222222222222222*x1)/(pow(x0,2) + pow(x1,2)) +
    0.3950617283950617*arg(9. + 16.*(x0 + complex<double>(0,1)*x1)) - 0.3950617283950617*arg(x0 + complex<double>(0,1)*x1);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}
int alpha_df(const gsl_vector *x, void *params, gsl_matrix *J) {
  double a = ((rparams *) params)->a;
  double b = ((rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);

  const double df00 = 0. + (0.2222222222222222*pow(x0,2))/pow(pow(x0,2) + pow(x1,2),2) - 
    (0.2222222222222222*pow(x1,2))/pow(pow(x0,2) + pow(x1,2),2) - 
    (0.3950617283950617*x0)/(pow(x0,2) + pow(x1,2)) + 56.888888888888886/(pow(9. + 16.*x0,2) + 256.*pow(x1,2)) + 
    (101.1358024691358*x0)/(pow(9. + 16.*x0,2) + 256.*pow(x1,2));
  const double df01 = (-0.4444444444444444*x0*x1)/pow(pow(x0,2) + pow(x1,2),2) + (0.3950617283950617*x1)/(pow(x0,2) + pow(x1,2)) - 
    (101.1358024691358*x1)/(pow(9. + 16.*x0,2) + 256.*pow(x1,2));
  const double df10 = -df01;
  const double df11 = df00;

  gsl_matrix_set (J, 0, 0, df00);
  gsl_matrix_set (J, 0, 1, df01);
  gsl_matrix_set (J, 1, 0, df10);
  gsl_matrix_set (J, 1, 1, df11);

  return GSL_SUCCESS;
}
int alpha_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
  alpha_f(x, params, f);
  alpha_df(x, params, J);

  return GSL_SUCCESS;
}

int rosenbrock_f(const gsl_vector *x, void *params, gsl_vector *f) {
  double a = ((rparams *) params)->a;
  double b = ((rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);

  const double y0 = a*(1 - x0);
  const double y1 = b*(x1 - pow(x0, 2));

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}
int rosenbrock_df(const gsl_vector *x, void *params, gsl_matrix *J) {
  double a = ((rparams *) params)->a;
  double b = ((rparams *) params)->b;

  const double x0 = gsl_vector_get (x, 0);

  const double df00 = -a;
  const double df01 = 0;
  const double df10 = -2*b*x0;
  const double df11 = b;

  gsl_matrix_set (J, 0, 0, df00);
  gsl_matrix_set (J, 0, 1, df01);
  gsl_matrix_set (J, 1, 0, df10);
  gsl_matrix_set (J, 1, 1, df11);

  return GSL_SUCCESS;
}
int rosenbrock_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
  rosenbrock_f(x, params, f);
  rosenbrock_df(x, params, J);

  return GSL_SUCCESS;
}

void print_state(size_t iter, gsl_multiroot_fdfsolver *s) {
  printf ("iter = %3u x = % .10f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

int main (int argc, char* argv[]) {
  // test multidimensional roots-finding
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 2;
  struct rparams p = {1.0, 10.0};
  gsl_multiroot_function_fdf f = {&alpha_f, &alpha_df, &alpha_fdf, n, &p};

  double x_init[2] = { 0.07, -0.02 };
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, x_init[0]);
  gsl_vector_set(x, 1, x_init[1]);

  T = gsl_multiroot_fdfsolver_hybridsj;
  s = gsl_multiroot_fdfsolver_alloc(T, n);
  gsl_multiroot_fdfsolver_set(s, &f, x);

  print_state(iter, s);

  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s);
    print_state(iter, s);

    if(status)
      break;

    status = gsl_multiroot_test_residual (s->f, 1e-7);
  } while( status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror(status));

  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);


  // cout.precision(17);
  // string outputFilePath = "./output/fits.csv";
  // if (argc > 1) {
  //   printf("Error no outputfile argument. Try ./build/FESR ./output/folder/fits.dat");
  //   outputFilePath = argv[1];
  // }

  // std::ifstream configFile("./configuration.json");
  // json config;
  // configFile >> config;

  // const Constants constants(config);
  // const Chisquared chisquared(config, constants);

  // // Numerics num(constants);
  // // cout << num.complexContourIntegral(testFunction) << endl;

  // // AdlerFunction adler(4, constants);
  // // cout << adler.D0(3.0, 3.0, 0.32307, 5) << endl;
  // // cout << 2.0*adler.D0CInt(3.0, Weight(config["parameters"]["weight"].get<int>()), 0.32307 , 5) << endl;
  // // cout << 2.0*adler.D0CInt(3.1572314596, Weight(config["parameters"]["weight"].get<int>()), 0.32307 , 5) << endl;

  // AlphaS amu(constants, 5);
  // complex<double> s(3.1355190333342473,-0.3696361468885539);
  // complex<double> mu2(-3.1355190333342473,0.3696361468885539);
  // cout << amu(mu2, constants.kSTau, 0.3179/constants.kPi) << endl;



  // cout << "th" << endl;
  // TheoreticalMoments thMom(config);
  // cout << thMom.cIntVpAD0FO(3.1572314596, Weight(1), 0.32307, 5) << endl;

  // cout << "Are ExpMom same: " << compareVectors(readExpMomentsFromFile(), chisquared.expMom_.getExpPlusPionPoleMoments()) << endl;
  // cout << "Is CovMat same: " << compareSquareMatrix(readMatrixFromFile(9, "./data/covMat.dat"), chisquared.expMom_.covarianceMatrix) << endl;
  // Numerics num(constants);

  // cout << "invInvCovMat: " << prod(chisquared.expMom_.covarianceMatrix, chisquared.expMom_.inverseCovarianceMatrix) << endl;


  // beta 4th order
  // cout << "Matthias CHI2: " << chisquared(0.32307175541329564, 2.1e-2, -0.309486307083497, -3.0869411117237483e-2) << endl;
  // cout << "MY CHI2: " << chisquared(0.32306509503077424, 2.1e-2, -0.30906464837998099, -0.029781360856290042)  << endl;

  // beta 5th order
  // cout << "Matthias CHI2: " << chisquared(0.32326096168471358, 2.1e-2, -0.31488720134123538, -2.6524803026353995e-2) << endl;
  // cout << "MY CHI2: " << chisquared(0.31541087265066009, 2.1e-2, 0.81120226780379001, 5.2853783566688879)  << endl;


  // compare with matthias
  // cout << "chi2Mat \t" << chisquared(0.32136770578073276, 2.1e-2, -0.30949, -3.0869e-2) << endl;
  // chisquared.log(0.32136770578073276, 2.1e-2, -0.30949, -3.0869e-2);


  // MINUIT
  // Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // // set tolerances
  // // min->SetMaxFunctionCalls(10000000); // for Minuit2
  // // min->SetMaxIterations(10000000); // for GSL
  // min->SetTolerance(1.0);

  // min->SetStrategy(2);
  // min->SetPrintLevel(3); // activate logging

  // // function wrapper
  // Functor chi2(chisquared, 4);

  // min->SetFunction(chi2);

  // // set free variables to be minimized
  // if (config["variables"]["astau"]["fixed"]) {
  //   min->SetFixedVariable(0, "astau", config["variables"]["astau"]["value"]);
  // } else {
  //   min->SetVariable(0, "astau", config["variables"]["astau"]["value"], config["variables"]["astau"]["stepSize"]);
  // }
  // if (config["variables"]["aGGInv"]["fixed"]) {
  //   min->SetFixedVariable(1, "aGGInv", config["variables"]["aGGInv"]["value"]);
  // } else {
  //   min->SetVariable(1, "aGGInv", config["variables"]["aGGInv"]["value"], config["variables"]["aGGInv"]["stepSize"]);
  // }
  // if (config["variables"]["rhoVpA"]["fixed"]) {
  //   cout << "fixed" << endl;
  //   min->SetFixedVariable(2, "rhoVpA", config["variables"]["rhoVpA"]["value"]);
  // } else {
  //   min->SetVariable(2, "rhoVpA", config["variables"]["rhoVpA"]["value"], config["variables"]["rhoVpA"]["stepSize"]);
  // }
  // if (config["variables"]["c8VpA"]["fixed"]) {
  //   min->SetFixedVariable(3, "c8VpA", config["variables"]["c8VpA"]["value"]);
  // } else {
  //   min->SetVariable(3, "c8VpA", config["variables"]["c8VpA"]["value"], config["variables"]["c8VpA"]["stepSize"]);
  // }

  // // minimize!
  // min->Minimize();
  // const double *xs = min->X();
  // const double *errors = min->Errors();
  // const double edm = min->Edm();
  // chisquared.log(xs[0], xs[1], xs[2], xs[3]);
  // const double chi2AtMin = chisquared(xs[0], xs[1], xs[2], xs[3]);
  // // min->PrintResults();

  // writeOutput(outputFilePath, xs, errors, chi2AtMin, edm, config);

  return 0;
}
