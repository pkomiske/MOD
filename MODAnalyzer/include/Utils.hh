#ifndef __UMOD_UTILS__
#define __UMOD_UTILS__

// standard library
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>

// external library headers
#include "fastjet/PseudoJet.hh"

// internal headers
#include "uMODAnalyzer.hh" // for mod::Particle

// mathematical constants
const double PI    = 3.141592653589793238462643383279502884197;
const double TWOPI = 2*PI;

// basic typedefs
typedef unsigned int uint;
typedef std::vector<int> Vint;
typedef std::vector<uint> Vuint;

// fastjet typedefs
typedef std::vector<fastjet::PseudoJet> VPJ;

// uMODAnalyzer typedefs
typedef std::vector<mod::Particle> VPart;
typedef std::vector<mod::AK5> VAK5;

// make particle phis close to jet phi
inline double phi_fix(double phi, double ref_phi) {
  double diff = phi - ref_phi;
  if (diff > PI) phi -= TWOPI;
  else if (diff < -PI) phi += TWOPI;
  return phi; 
}

// round masses
inline std::string m_round(double m) {
  std::ostringstream ss;

  // note that these values are somewhat carefully tuned from pythia masses
  if (m < 1e-5) m = 0;
  ss << std::setprecision(5) << m;

  return ss.str();
}

// calculate rapidity
inline double rap_from_PtEtaM(double pt, double eta, double m) {
  double m2(m*m), pt2(pt*pt), se(sinh(eta)), ce(cosh(eta));
  return log((sqrt(m2 + pt2*ce*ce) + pt*se)/(sqrt(m2 + pt2)));
}

// more kinematic functions
fastjet::PseudoJet apply_jec_fixed_eta(const fastjet::PseudoJet & pj, double jec);
VPJ apply_jecs(const VPJ & pjs, const std::vector<double> & jecs);

// pseudojet creation from mod structs
VPJ to_pjs(const VPart & particles);
VPJ to_pjs(const VAK5 & ak5s);

// from Table V of https://arxiv.org/pdf/1704.05842.pdf
uint jet_quality(const double neutralHadronFraction,
                 const double neutralEMFraction,
                 const double chargedHadronFraction,
                 const double chargedEMFraction,
                 const uint   nConsts,
                 const uint   nChargedConsts,
                 const double eta);

// argsorting jets by pt
Vuint argsort_by_pt(const VPJ pjs); 
void argsort(Vuint & indices, const std::vector<double> & values);

class IndexedSortHelper {
public:

  inline IndexedSortHelper (const std::vector<double> & values) : values_(values) {}
  inline int operator() (const int i1, const int i2) const {
    return  values_[i1] < values_[i2];
  }

private:
  const std::vector<double> & values_;
};

#endif // __UMOD_UTILS__
