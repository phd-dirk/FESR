#include "utils.hpp"

void addCSV(double value, ofstream &outputFile, bool last = false) {
  if(!last) {
    outputFile << value << ",";
  } else {
    outputFile << value;
  }
}

string variablesOutput(Minimizer* min) {
  const double *xs = min->X();
  const double *errors = min->Errors();
  string out = "";

  for(uint i = 0; i < min->NDim(); i++) {
    const bool isFixed = min->IsFixedVariable(i);
    // name
    out += min->VariableName(i) + "\t=\t" + std::to_string(xs[i]);

    // error
    if(!isFixed) {
      out += "\t+/-\t" + std::to_string(errors[i]) + "\t";
    }

    // fixed
    if(isFixed) {
      out += "\t(fixed)";
    }

    out += "\n";
  }
  out += "\n";
  return out;
}

string correlationOutput(Minimizer* min) {
  string out;
  out += "Correlation Matrix: \n";
  for(uint i = 0; i < min->NDim(); i++) {
    for(uint j = 0; j < min->NDim(); j++) {
      out += std::to_string(min->Correlation(i, j));
      if(j < min->NDim() - 1) out += ", ";
    }
    out += "\n";
  }
  out += "\n";
  return out;
}

string deltaOutput(Minimizer* min, const Configuration &config) {
  string out;
  ThMoms thMom(config);
  const double *xs = min->X();

  out += "Deltas (FO): \n";
  out += "delta_V+A^(0) \t" +
      std::to_string(
          thMom.del0(
              config.sTau_, config.inputs_[0].weight, config.sTau_, xs[0],
              Configuration::adlerCoefficients(
                3, Configuration::betaCoefficients(3, 3)
              ),
              config.order_
          )
      ) + "\n";
  out += "delta_V+A^(4) \t" +
      std::to_string(
          thMom.del4(
              config.sTau_, config.inputs_[0].weight, config.sTau_, xs[0], config.aGGInv.value
          )
      ) + "\n";
  out += "delta_V+A^(6) \t" +
      std::to_string(
          thMom.del6(
              config.sTau_, config.inputs_[0].weight, config.rhoVpA.value
          )
      ) + "\n";
  out += "delta_V+A^(8) \t" +
      std::to_string(
          thMom.del8(
              config.sTau_, config.inputs_[0].weight, config.c8VpA.value
          )
      ) + "\n";
  out += "delta_V+A^(10) \t" +
    std::to_string(
      thMom.del10(
        config.sTau_, config.inputs_[0].weight, config.c10.value
      )
    ) + "\n";
  out += "delta_V+A^(12) \t" +
    std::to_string(
      thMom.del12(
        config.sTau_, config.inputs_[0].weight, config.c12.value
      )
    ) + "\n";
  out += "\n";
  return out;
}

void writeToCSV(
    string &csvFilePath,
    Minimizer* min,
    const Configuration config
) {
  ofstream csvFile;
  csvFile.open(csvFilePath, std::ios::app);
  csvFile << std::setprecision(15);

  if(csvFile.is_open())
  {
    const double *xs = min->X();
    const double *errors = min->Errors();
    // min status
    addCSV(min->Status(), csvFile);
    // alpha_s
    addCSV(xs[0], csvFile);
    addCSV(errors[0], csvFile);
    // aGGInv
    addCSV(xs[1], csvFile);
    addCSV(errors[1], csvFile);
    // O6
    addCSV(xs[2], csvFile);
    addCSV(errors[2], csvFile);
    // O8
    addCSV(xs[3], csvFile);
    addCSV(errors[3], csvFile);

    // write FVAL
    addCSV(min->MinValue(), csvFile);
    // dof
    addCSV(config.dof(), csvFile);
    addCSV(min->MinValue()/config.dof(), csvFile);
    addCSV(min->Edm(), csvFile, true);

    csvFile << endl;
    csvFile.close();
  }
}

void writeToOutputFile(
    const string &configFilePath,
    Minimizer* min,
    const Configuration &config
) {
  // std::size_t pos = configFilePath.find_last_of("/\\");
  std::size_t pos = configFilePath.find_last_of(".");
  std::string outputFilePath = configFilePath.substr(0,pos) + ".txt";

  ofstream outputFile;
  outputFile.open(outputFilePath);
  outputFile << std::setprecision(15);

  ifstream configFile;
  configFile.open(configFilePath);

  if(outputFile.is_open() && configFile.is_open())
  {
    outputFile << "Minuit2Minimizer\t Status:\t" << min->Status() << endl;
    outputFile << "FVAL\t=\t" << min->MinValue() << endl;
    outputFile << "Edm\t=\t" << min->Edm() << endl;
    outputFile << "dof = \t" << config.dof() << std::endl;
    outputFile << "chi2/dof = \t" << min->MinValue()/config.dof() << std::endl;
    outputFile << "Nfcn\t=\t" << min->NCalls() << endl;

    const double *xs = min->X();
    const double *errors = min->Errors();

    // variables
    outputFile << variablesOutput(min);

    // correlation matrix
    outputFile << correlationOutput(min);

    // deltas
    outputFile << deltaOutput(min, config);


    // add configFile
    outputFile << configFile.rdbuf();

    configFile.close();
    outputFile.close();
  }
}

void writeOutput (
    const string configFilePath,
    Minimizer* min,
    const Configuration config
) {
  std::size_t pos = configFilePath.find_last_of("/\\");
  std::string csvFilePath = configFilePath.substr(0,pos) + "/fits.csv";

  writeToCSV(csvFilePath, min, config);
  writeToOutputFile(configFilePath, min, config);
}


void writeOutput(const string text, const string outputFilePath)
{
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::app);
  outputFile << text << endl;
}

int utils::momCount(std::vector<Input> inputs)
{
  int momCount = 0;
  for(auto const &input: inputs) {
    momCount += input.s0s.size();
  }
  return momCount;
}
