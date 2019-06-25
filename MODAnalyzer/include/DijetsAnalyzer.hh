#ifndef __UMOD_JETS_ANALYZER__
#define __UMOD_JETS_ANALYZER__

// standard library
#include <string>
#include <vector>

// external libraries
#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"

// internal headers
#include "uMODAnalyzer.hh"
#include "JetMatcher.hh"
#include "Utils.hh"

// constants
#define N_DECAY_ATTEMPTS 15

class DijetsAnalyzer : public mod::Analyzer {
public:

  // constructor with fully explicit arguments
  DijetsAnalyzer(std::string & inFilePath, std::string & outFilePath, std::string & fileType,
                 int num,                  int printEvery,            int outprecision,
                 bool doDecays,            bool trimTrigNames,        bool countTriggers,
                 double jet_R,             double minJetPt,           double minJetPt4PFCs,
                 double matchCMSPtFrac,    double matchCMSDeltaR,
                 double matchGenPtFrac,    double matchGenDeltaR,
                 double matchHardPtFrac,   double matchHardDeltaR);

  // pythia related functions
  void setup_pythia();
  VPart decay_particles(const VPart & particles, const uint nAttempts = N_DECAY_ATTEMPTS) const;

  // fastjet related functions
  VPJ cluster_pfcs(const VPart & pfcs);
  VPJ cluster_gens(const VPart & gens);
  VPJ cluster_raw_ak5_jets(const VPJ & pjs, fastjet::ClusterSequence * & cs, double minJetPt) const;

  // processing related functions
  std::vector<bool> get_output_mask(uint n, const Vint & ref2new) const;
  void process_event(const mod::Event & event);
  void print_progress(std::ostream & ostr);

  // output related functions
  void output_event_info(const mod::Event & event);
  void output_ak5s(const VPJ & pjs, const VAK5 & ak5s);

  void output_jets_masked(const VPJ & jets, const std::vector<bool> & mask, 
                          const VPart & particles, const char type);

  void output_particle(const fastjet::PseudoJet & pj, const mod::Particle & p, 
                       const char type, const int i, const double phi_ref);

  void output_match_inds(const Vint & ref2new, const char * names); 
  void output_lumiblocks();

private:

  // variables for redecaying
  const bool doDecays_, trimTrigNames_, countTriggers_;
  Pythia8::Pythia * pythia_;

  // variables related to jet finding
  const double jet_R_, minJetPt_, minJetPt4PFCs_;
  fastjet::JetDefinition jetDef_;
  fastjet::ClusterSequence * csPFC_, * csGen_;
  VPJ pfcPJs_, genPJs_;

  // objects for matching various pairs of jets
  int verbosity_;
  JetMatcher cmsMatcher_, genMatcher_, hardMatcher_;

  // counters
  uint numNoCMSJets_, numPFCsMissing_, numPythiaFail_;

};

#endif // __UMOD_JETS_ANALYZER__
