#include "EventKaons.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int ac, char *av[])
{
  // Declare the CLI options and parameters.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("input-file", po::value<string>(), "Path to input ROOT nuptle")
      ("help,h", "produce this help message and exit")
      ("name,n", po::value<string>()->default_value("unnamed"),
          "Name of the processed data");
  
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
  
  // Process ntuple: generate plots and dataset.
  EventKaons kaon_distribution(
      vm["input-file"].as<string>(),
      vm["name"].as<string>(),
      "|p_{0}| + |p_{1}|", "GeV/c", 2.0, 5.25);

  // Process raw data.
  kaon_distribution.processNtuple();

  // Save processed data.
  kaon_distribution.drawHistograms();
  kaon_distribution.saveHistograms();
  
  return 0;
}
