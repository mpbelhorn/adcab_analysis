#include "Fit3D.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int ac, char *av[])
{
  bool generate_flag = false;
  bool fit_flag = false;
  
  // Declare the CLI options and parameters.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("input-file", po::value<string>(), "Path to input ROOT nuptle")
      ("help,h", "produce this help message and exit")
      ("name,n", po::value<string>()->default_value("unnamed"),
          "Name of the processed data")
      ("generate,g", po::bool_switch(&generate_flag),
          "Generate PDFs from input data")
      ("fit,f", po::bool_switch(&fit_flag), "Fit data to generated PDFs");
  
  po::positional_options_description pos;
  pos.add("input-file", 1);
  
  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).
      options(desc).positional(pos).run(), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << "Usage: Fit3D [options] input-file\n";
    std::cout << desc << "\n";
    return 1;
  }

  if (vm.count("input-file") != 1) {
    std::cout << "Input file specified incorrectly!\n";
    std::cout << "Usage: Fit3D [options] input-file\n";

    return 1;
  }
  
  if (!(fit_flag || generate_flag)) {
    std::cout << "Specify at least one of the flags 'fit' or 'generate'.\n";
    return 1;
  }
  
  // Process ntuple: generate plots and dataset.
  Fit3D fitter(
      vm["input-file"].as<string>(),
      vm["name"].as<string>(),
      "|p_{0}| + |p_{1}|", "GeV/c", 2.0, 5.25,
      "Cos(#theta_{ll})", "", -0.80, 0.98,
      "#Delta z", "cm", 0.0, 0.18);

  // Process raw data.
  fitter.setCreateDataSet(true);
  fitter.processNtuple();

  // Save processed data.
  fitter.drawHistograms();
  fitter.saveHistograms();
  fitter.saveDataSet();
  
  // Generate models, if requested.
  if (generate_flag) {
    fitter.generateModels(2);
  }
  
  // Fit data, if requested.
  if (fit_flag) {
    fitter.fitData();
    fitter.plotAsymmetry();
  }
  
  return 0;
}
