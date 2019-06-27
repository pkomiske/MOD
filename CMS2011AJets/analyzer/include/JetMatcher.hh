#ifndef __UMOD_JET_MATCHER__
#define __UMOD_JET_MATCHER__

// standard library
#include <sstream>
#include <string>

// external library
#include "fastjet/PseudoJet.hh"

// internal headers
#include "Utils.hh"

class JetMatcher {
public:

  // constructor
  JetMatcher(double ptFrac_, double deltaRMax_, 
             std::string newName_, std::string refName_, std::string filename_, 
             int verbosity);

  ~JetMatcher() {}

  // function that actual determines if two jets match
  inline bool match_specific_jets(const fastjet::PseudoJet & jet1, const fastjet::PseudoJet & jet2) const {
    
    // kinematic considerations
    if (jet1.delta_R(jet2) > deltaRMax || 2*abs(jet1.pt() - jet2.pt())/(jet1.pt() + jet2.pt()) > ptFrac) 
      return false;
    
    return true;
  }

  // match collections of jets
  Vint match(const VPJ & newJets, const VPJ & refJets, uint event_i);

  // output counters
  std::string stats() const;
  inline std::string cols() const { return cols_; }

  // convenience
  Vint operator()(const VPJ & newJets, const VPJ & refJets, uint event_i) {
    return match(newJets, refJets, event_i);
  }

  inline double pt_frac() const { return ptFrac; }
  inline double delta_r_max() const {return deltaRMax; }

private:

  // whoami
  const std::string newName, refName, filename;

  // matching parameters
  const double ptFrac, deltaRMax;

  // verbosity
  const int verbosity_;

  // counters
  uint numNewMatchedMultipleTimes, numRefMatchedMultipleTimes,
       numHardestTwoNewFailMatch, numHardestTwoRefFailMatch,
       numHardestTwoRefMatchedNotHardestTwoNew;

  // stat columns
  const std::string cols_;

};

#endif // __UMOD_JET_MATCHER__
