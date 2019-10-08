#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "MODDataset.hh"

#define OUTPRECISION 14

int main(int argc, char** argv) {

  // command line arguments
  if (argc < 4) {
    std::cerr << "Usage: ./MODDatasetTest MOD_HDF5_FILENAME OUTPATH INDEX\n";
    return 1;
  }

  // initialize dataset
  MOD::Dataset dataset(argv[1]);

  // initialize output text files
  std::string outpath(argv[2]), index(argv[3]);
  std::ofstream of_jets_i(outpath + "/jets_i_" + index + ".txt"),
                of_jets_f(outpath + "/jets_f_" + index + ".txt"),
                of_pfcs(outpath + "/pfcs_" + index + ".txt"),
                of_gens(outpath + "/gens_" + index + ".txt");

  // set appropriate precision for floats
  of_jets_f << std::setprecision(OUTPRECISION);
  of_pfcs << std::setprecision(OUTPRECISION);
  of_gens << std::setprecision(OUTPRECISION);

  // iterate through events
  while (dataset.next()) {

    // output jets_i
    for (auto i : dataset.jets_i()) of_jets_i << i << ' ';
    of_jets_i << '\n';

    // output jets_f
    for (auto f : dataset.jets_f()) of_jets_f << f << ' ';
    of_jets_f << '\n';

    // output pfcs if we have them
    if (dataset.hasPFCs()) {
      for (auto & p : dataset.pfcs())
        of_pfcs << p.pt     << ' '
                << p.y      << ' '
                << p.phi    << ' '
                << p.m      << ' '
                << p.pid    << ' '
                << p.vertex << '\n';
      of_pfcs << '\n';
    }

    // output gens if we have them
    if (dataset.hasGENs()) {
      for (auto & g : dataset.gens())
        of_gens << g.pt     << ' '
                << g.y      << ' '
                << g.phi    << ' '
                << g.m      << ' '
                << g.pid    << ' '
                << g.vertex << '\n';
      of_gens << '\n';
    }
  }

  return 0;
}