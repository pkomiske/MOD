// standard library
#include <cstdlib>
#include <iostream>
#include <string>

// internal headers
#include "CmdLine.hh"
#include "DijetsAnalyzer.hh"

int main(int argc, char** argv) {

  // parse command line
  CmdLine cmdline(argc, argv);

  // get required positional arguments
  if (argc < 4) {
    std::cerr << "Usage: dijets_analyzer INFILEPATH OUTFILEPATH FILETYPE [OPTIONS]\n";
    return 1;
  }

  std::string inFilePath(argv[1]), outFilePath(argv[2]), fileType(argv[3]);

  // integer arguments
  int num          = cmdline.value<int>("-n",            -1),
      printEvery   = cmdline.value<int>("-print",        -1),
      outprecision = cmdline.value<int>("-outprecision", 10);

  // boolean arguments
  bool doDecays      = !cmdline.present("-no-decays"),
       trimTrigNames = !cmdline.present("-keeptrignames"),
       countTriggers = !cmdline.present("-no-trigcounts");

  // jet parameters
  double jet_R         = cmdline.value<double>("-jetR",          0.5),
         minJetPt      = cmdline.value<double>("-minjetpt",      2.5),
         minJetPt4PFCs = cmdline.value<double>("-minjetpt-pfcs", 300);

  // matching parameters
  double matchCMSPtFrac  = cmdline.value<double>("-match-cms-ptfrac",       0.000001),
         matchCMSDeltaR  = cmdline.value<double>("-match-cms-deltar",       0.000001),
         matchGenPtFrac  = cmdline.value<double>("-match-gen-ptfrac",       1000),
         matchGenDeltaR  = cmdline.value<double>("-match-gen-deltar-mult",  1) * jet_R,
         matchHardPtFrac = cmdline.value<double>("-match-hard-ptfrac",      1000),
         matchHardDeltaR = cmdline.value<double>("-match-hard-deltar-mult", 2) * jet_R;

  // ensure no unexpected options
  cmdline.assert_all_options_used(); 

  // run analyzer
  DijetsAnalyzer(inFilePath,      outFilePath,    fileType,
                 num,             printEvery,     outprecision,
                 doDecays,        trimTrigNames,  countTriggers,
                 jet_R,           minJetPt,       minJetPt4PFCs,
                 matchCMSPtFrac,  matchCMSDeltaR,
                 matchGenPtFrac,  matchGenDeltaR,
                 matchHardPtFrac, matchHardDeltaR);

  return 0;
}
