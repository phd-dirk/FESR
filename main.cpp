#include <chrono>
#include <iostream>
#include <stdexcept>
#include "./src/configuration.hpp"
#include "./src/minuit.hpp"
#include "./src/utils.hpp"

#include "./src/chisquared.hpp"
#include <vector>

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


    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6  <<std::endl;

    return 0;
  }
  catch (const std::exception& e) {
    cout << e.what() << endl;
    return 1;
  }
}
