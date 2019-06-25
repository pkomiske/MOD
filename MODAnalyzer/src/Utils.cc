#include <cmath>
#include <sstream>

#include "Utils.hh"

fastjet::PseudoJet apply_jec_fixed_eta(const fastjet::PseudoJet & pj, double jec) {

  fastjet::PseudoJet newPJ(pj);
  double pt(pj.pt() * jec), eta(pj.eta()), phi(pj.phi()), m(pj.m());

  newPJ.reset_momentum_PtYPhiM(pt, rap_from_PtEtaM(pt, eta, m), phi, m);
  return newPJ;
}

VPJ apply_jecs(const VPJ & pjs, const std::vector<double> & jecs) {

  VPJ pjsCorrected;
  for (uint i = 0; i < pjs.size(); i++)
    pjsCorrected.push_back(apply_jec_fixed_eta(pjs[i], jecs[i]));

  return pjsCorrected;
}

VPJ to_pjs(const VPart & particles) {

  VPJ pjs(particles.size());
  for (uint i = 0; i < particles.size(); i++) {
    const mod::Particle & particle(particles[i]);
    pjs[i].reset_momentum(particle.px, particle.py, particle.pz, particle.e);
    pjs[i].set_user_index(i);
  }

  return pjs;
}

VPJ to_pjs(const VAK5 & ak5s) {

  VPJ pjs(ak5s.size());
  for (uint i = 0; i < ak5s.size(); i++) {
    const mod::AK5 & ak5(ak5s[i]);
    double y(rap_from_PtEtaM(ak5.pt, ak5.eta, ak5.m));

    pjs[i].reset_momentum_PtYPhiM(ak5.pt, y, ak5.phi, ak5.m);
    pjs[i].set_user_index(i);
  }

  return pjs;
}

// from Table V of https://arxiv.org/pdf/1704.05842.pdf
uint jet_quality(const double neutralHadronFraction,
                 const double neutralEMFraction,
                 const double chargedHadronFraction,
                 const double chargedEMFraction,
                 const uint   nConsts,
                 const uint   nChargedConsts,
                 const double eta) {

  if (fabs(eta) < 2.4) {
    if (neutralHadronFraction < 0.90 &&
        neutralEMFraction     < 0.90 &&
        nConsts               > 1    &&
        chargedHadronFraction > 0    &&
        chargedEMFraction     < 0.99 &&
        nChargedConsts        > 0) return 3;

    else if (neutralHadronFraction < 0.95 &&
             neutralEMFraction     < 0.95 &&
             nConsts               > 1    &&
             chargedHadronFraction > 0    &&
             chargedEMFraction     < 0.99 &&
             nChargedConsts        > 0) return 2;

    else if (neutralHadronFraction < 0.99 &&
             neutralEMFraction     < 0.99 &&
             nConsts               > 1    &&
             chargedHadronFraction > 0    &&
             chargedEMFraction     < 0.99 &&
             nChargedConsts        > 0) return 1;
  }

  else {
    if (neutralHadronFraction < 0.90 &&
        neutralEMFraction     < 0.90 &&
        nConsts               > 1) return 3;

    else if (neutralHadronFraction < 0.95 &&
             neutralEMFraction     < 0.95 &&
             nConsts               > 1) return 2;

    else if (neutralHadronFraction < 0.99 &&
             neutralEMFraction     < 0.99 &&
             nConsts               > 1) return 1;
  }

  return 0;

  /*// 2 indicates always bad
  if (nConsts <= 1 || // bad because too few constituents
      (fabs(eta) < 2.4 && // this jet must pass charged requirements also
        (chargedHadronFraction == 0. || chargedEMFraction >= 0.99 || nChargedConsts == 0))) 
    return 2;

  // return max value of neutral fractions (this always determines quality at this point)
  return std::max(neutralHadronFraction, neutralEMFraction);*/
}

void argsort(Vuint & indices, const std::vector<double> & values) {
  IndexedSortHelper indexedSortHelper(values);
  std::sort(indices.begin(), indices.end(), indexedSortHelper);
}

Vuint argsort_by_pt(const VPJ pjs) {

  uint numPJs(pjs.size());

  Vuint indices(numPJs);
  std::vector<double> negPt2s(numPJs);
  for (uint i = 0; i < numPJs; i++) {
    indices[i] = i;
    negPt2s[i] = -pjs[i].kt2();
  }

  argsort(indices, negPt2s);

  return indices;
}
