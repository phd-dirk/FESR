#include <chrono>
#include <iostream>
#include <stdexcept>
#include "./src/configuration.hpp"
#include "./src/minuit.hpp"
#include "./src/utils.hpp"

#include "./src/chisquared.hpp"
#include <vector>

// #include "./src/experimentalMoments.hpp"
// #include <boost/numeric/ublas/matrix.hpp>

int main (int argc, char* argv[]) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::cout.precision(17);

  try {
    string configFilePath = "./configuration.json";
    if (argc == 2) {
      configFilePath = argv[1];
    }
    const Configuration config(configFilePath);


    // SpecEnd
    // const int num_bins = 6;
    // const ROOT::Math::Minimizer* min = Minuit::spec_end(config, num_bins);

    // FESR
    const ROOT::Math::Minimizer* min = Minuit::FESR(config);
    writeOutput(configFilePath, min, config);


    // TEST
    // ExpMoms expMom(
    //   "/Users/knowledge/Developer/PhD/FESR/aleph.json",
    //   {
    //     {
    //       Weight(14),
    //       {
    //         3.1572314596000002, 3.0, 2.800, 2.600, 2.4, 2.3, 2.2
    //       }
    //     }
    //   },
    //   3.1572314596000002, // sTau
    //   17.815, // be
    //   0.023, // dBe
    //   0.97420, // Vud
    //   0.00021, // dVud
    //   1.0198, // SEW
    //   0.0006, // dSEW
    //   92.21e-3, // fPi
    //   0.14e-3, // dFPi
    //   0.13957018, // pionMinusMass
    //   0.99743669 // RVANormalization
    // );

    // std::vector<double> moments = expMom();
    // for(auto &mom: moments) {
    //   std::cout << mom << std::endl;
    // }

    // // boost::numeric::ublas::matrix<double> C;
    // // C = boost::numeric::ublas::prod(expMom.covMat_, expMom.invCovMat_);

    // std::cout << expMom.invCovMat_ << std::endl;


    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6  <<std::endl;

    return 0;
  }
  catch (const std::exception& e) {
    cout << e.what() << endl;
    return 1;
  }
}
