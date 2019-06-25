#include "JetMatcher.hh"

JetMatcher::JetMatcher(double ptFrac_, double deltaRMax_,
                       std::string newName_, std::string refName_, std::string filename_,
                       int verbosity) :
// whoami
newName(newName_), refName(refName_), filename(filename_),

// matching parameters
ptFrac(ptFrac_), deltaRMax(deltaRMax_), 

// printing
verbosity_(verbosity),

// counters
numNewMatchedMultipleTimes(0), numRefMatchedMultipleTimes(0),
numHardestTwoNewFailMatch(0), numHardestTwoRefFailMatch(0),
numHardestTwoRefMatchedNotHardestTwoNew(0),
cols_("          h2RefNoMatch   h2NewNoMatch   h2Ref!~h2New  nRefMatchMany  nNewMatchMany")
{}

// ref should be length at most 2 and new can be longer
Vint JetMatcher::match(const VPJ & newJets, const VPJ & refJets, const uint event_i) {

  uint sizeNew(newJets.size()), sizeRef(refJets.size());

  // containers of indices with the mapping between new and old jets
  Vint new2ref(sizeNew, -1), ref2new(sizeRef, -1);

  for (uint j = 0; j < sizeRef; j++) {
    const fastjet::PseudoJet & refJet(refJets[j]);

    for (uint i = 0; i < sizeNew; i++) {
      const fastjet::PseudoJet & newJet(newJets[i]);

      // if new matches ref and ref has not previously matched
      if (match_specific_jets(newJet, refJet)) {
        if (ref2new[j] == -1) {
          if (new2ref[i] == -1) {
            new2ref[i] = j;
            ref2new[j] = i;
            break;
          }
          else {
            numNewMatchedMultipleTimes++;
            if (verbosity_ >= 1)
              std::cout << filename << ", event_i " << event_i << ": " << refName << " jet " << j << " with pt " << refJet.pt() 
                        << " GeV matched " << newName << " jet " << i << " with pt " << newJet.pt() 
                        << " that has already matched " << refName << " jet " << new2ref[i] << '\n';
          }
        }
        else {
          numRefMatchedMultipleTimes++;
          if (verbosity_ >= 1)
            std::cout << filename << ", event_i " << event_i << ": " << newName << " jet " << i << " with pt " << newJet.pt() 
                      << " GeV matched " << refName << " jet " << j << " with pt " << refJet.pt() 
                      << " that has already matched " << newName << " jet " << ref2new[j] << '\n';
        }
      }
    }
  }

  // check that highest two new jets matched
  for (uint i = 0, m = std::min<uint>(2, sizeNew); i < m; i++) {
    if (new2ref[i] == -1) {
      if (verbosity_ >= 1)
        std::cout << filename << ", event_i " << event_i << ": " << newName << " jet " << i 
                  << " did not match a " << refName << " jet!\n";
      numHardestTwoNewFailMatch++;
    }
  }

  // check that highest two ref jets matched and get matching new jets
  VPJ matchedNewJets;
  for (uint i = 0, m = std::min<uint>(2, sizeRef); i < m; i++) {
    int new_i(ref2new[i]);

    if (new_i == -1) {
      if (verbosity_ >= 1)
        std::cout << filename << ", event_i " << event_i << ": " << refName << " jet " << i 
                  << " did not match a " << newName << " jet!\n";
      numHardestTwoRefFailMatch++;
    }
    else if (new_i >= 2) {
      numHardestTwoRefMatchedNotHardestTwoNew++;
      if (verbosity_ >= 1)
        std::cout << filename << ", event_i " << event_i << ": " << refName << " jet " << i 
                  << " matched " << newName << " jet " << new_i << '\n';
    }
  }

  //return matchedNewJets;
  return ref2new;
}

std::string JetMatcher::stats() const {

  std::ostringstream oss;
  oss << newName << '-' << refName << std::setw(15) << numHardestTwoRefFailMatch
                                   << std::setw(15) << numHardestTwoNewFailMatch
                                   << std::setw(15) << numHardestTwoRefMatchedNotHardestTwoNew
                                   << std::setw(15) << numRefMatchedMultipleTimes
                                   << std::setw(15) << numNewMatchedMultipleTimes;
  return oss.str();
}
