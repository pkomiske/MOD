// standard library
#include <cmath>
#include <iomanip>
#include <iostream>

// internal headers
#include "DijetsAnalyzer.hh"

// this has to come after we import something from fastjet
#include "fastjet/version.hh"

DijetsAnalyzer::DijetsAnalyzer(std::string & inFilePath, std::string & outFilePath, std::string & fileType,
                               int num,                  int printEvery,            int outprecision,
                               bool doDecays,            bool trimTrigNames,        bool countTriggers,
                               double jet_R,             double minJetPt,           double minJetPt4PFCs,
                               double matchCMSPtFrac,    double matchCMSDeltaR,
                               double matchGenPtFrac,    double matchGenDeltaR,
                               double matchHardPtFrac,   double matchHardDeltaR) :

  // base class constructor
  Analyzer(inFilePath, outFilePath, fileType, num, printEvery, outprecision),

  // flags
  doDecays_(sim_ && doDecays),
  trimTrigNames_(trimTrigNames),
  countTriggers_(countTriggers),

  // fastjet
  jet_R_(jet_R), minJetPt_(minJetPt), minJetPt4PFCs_(minJetPt4PFCs),
  jetDef_(fastjet::antikt_algorithm, jet_R_),
  csPFC_(nullptr), csGen_(nullptr),

  // jet matching
  verbosity_(printEvery == -1 ? 0 : 1),
  cmsMatcher_(matchCMSPtFrac, matchCMSDeltaR, (sim_ ? "SIM" : "PFC"), "CMS", inFileName_, verbosity_),
  genMatcher_(matchGenPtFrac, matchGenDeltaR, "GEN", "SIM", inFileName_, verbosity_),
  hardMatcher_(matchHardPtFrac, matchHardDeltaR, "GEN", "HRD", inFileName_, verbosity_),

  // counters
  numNoCMSJets_(0), numPFCsMissing_(0), numPythiaFail_(0)
{
  assert(sim_ || fileType_ == "cms");

  if (doDecays_)
    setup_pythia();

  if (verbosity_ == 0)
    fastjet::ClusterSequence::set_fastjet_banner_stream(new std::ostringstream());

  parse_file();
  write_summary();
  output_lumiblocks();
  print_progress(outFile_);
}

void DijetsAnalyzer::setup_pythia() {

  uint seed;
  std::stringstream ss;
  ss << std::hex << ("0x" + inFileName_.substr(0, 8));
  ss >> seed;

  pythia_ = new Pythia8::Pythia("../share/Pythia8/xmldoc", verbosity_ > 0);

  pythia_->readString("ProcessLevel:all = off");
  pythia_->readString("HadronLevel:Decay = on");
  pythia_->readString("ParticleDecays:limitTau0 = on");
  pythia_->readString("ParticleDecays:tau0Max = 1000");
  pythia_->readString("Random:setSeed = on");
  pythia_->readString("Random:seed = " + std::to_string(seed % 900000000));
  pythia_->readString("Print:quiet = on");

  pythia_->init();
}

VPart DijetsAnalyzer::decay_particles(const VPart & particles, const uint nAttempts) const {

  VPart newParticles;
  uint n(0);
  do {

    // IMPORTANT: reset() adds a pseudoparticle to represent the event as a whole
    pythia_->event.reset();
      
    for (VPart::const_iterator it = particles.begin(), end = particles.end(); it != end; it++) {
      int pdgid(it->pdgid);
      double px(it->px), py(it->py), pz(it->pz), mass(pythia_->particleData.m0(pdgid));
      double e(sqrt(mass*mass + px*px + py*py + pz*pz));

      pythia_->event.append(pdgid, 91, 0, 0, px, py, pz, e, mass);
    }

    n++;
    if (n >= nAttempts) {
      std::cout << inFileName_ << ", event_i " << event_i_ 
                << ": Pythia decays failed after " << nAttempts << " attempts!\n";
      return newParticles;
    }

  } while (!pythia_->next());

  for (uint i = 0, size = pythia_->event.size(); i < size; i++) {
    Pythia8::Particle & p(pythia_->event[i]);
    
    uint idAbs = p.idAbs();
    if (p.isFinal() && idAbs != 12 && idAbs != 14 && idAbs != 16) 
      newParticles.emplace_back(p.px(), p.py(), p.pz(), p.e(), p.id());
  }

  return newParticles;
}

VPJ DijetsAnalyzer::cluster_pfcs(const VPart & pfcs) {
  pfcPJs_ = to_pjs(pfcs);
  return cluster_raw_ak5_jets(pfcPJs_, csPFC_, minJetPt_);
}

VPJ DijetsAnalyzer::cluster_gens(const VPart & gens) {
  genPJs_ = to_pjs(gens);
  return cluster_raw_ak5_jets(genPJs_, csGen_, minJetPt_);
}

inline VPJ DijetsAnalyzer::cluster_raw_ak5_jets(const VPJ & pjs, 
                                                fastjet::ClusterSequence * & cs, 
                                                double minJetPt) const {
  cs = new fastjet::ClusterSequence(pjs, jetDef_);
  return fastjet::sorted_by_pt(cs->inclusive_jets(minJetPt));
}

std::vector<bool> DijetsAnalyzer::get_output_mask(uint n, const Vint & ref2new) const {

  std::vector<bool> mask(n, false);

  for (uint i = 0; i < ref2new.size(); i++) {
    if (ref2new[i] != -1) 
      mask[ref2new[i]] = true;

    // ensure hardest two are output
    if (i < 2)
      mask[i] = true;
  }



  return mask;
}

void DijetsAnalyzer::process_event(const mod::Event & event) {

  output_event_info(event);

  // get pseudojet vectors of cms jets
  VPJ cmsAK5s(to_pjs(event.ak5s));

  // if no cms jets, continue
  if (cmsAK5s.size() == 0) {

    // record that we have no jets
    numNoCMSJets_++;
    return;
  }

  // get jecs of ak5 jets
  std::vector<double> jecs(event.ak5s.size());
  for (uint i = 0; i < jecs.size(); i++)
    jecs[i] = event.ak5s[i].jec;

  // get corrected cms ak5 jets
  VPJ cmsAK5sCorrected(apply_jecs(cmsAK5s, jecs));

  // find two hardest CMS jets, corrected and uncorrected
  Vuint cmsPtArgsort(argsort_by_pt(cmsAK5sCorrected));
  VPJ hardestTwoCMSJets, hardestTwoCMSJetsCorrected;
  for (uint i = 0; i < std::min<uint>(2, cmsPtArgsort.size()); i++) {
    hardestTwoCMSJets.push_back(cmsAK5s[cmsPtArgsort[i]]);
    hardestTwoCMSJetsCorrected.push_back(cmsAK5sCorrected[cmsPtArgsort[i]]);
  }

  // output hardest two AK5 jets
  output_ak5s(hardestTwoCMSJets, event.ak5s);

  // check that we have any particles at all
  if (event.pfcs.empty()) {

    // check that we have no jets with corrected pt over 300
    if (hardestTwoCMSJetsCorrected[0].pt() >= minJetPt4PFCs_) {
      std::cout << inFileName_ << ", event " << event_i_ 
                << ": has no particles but has a corrected jet over the threshold to write particles!\n";
      numPFCsMissing_++;
    }

    return;
  }

  // find jets among the pfcs
  VPJ pfcJets(cluster_pfcs(event.pfcs));

  // get matched pfc jets
  Vint cms2pfc(cmsMatcher_(pfcJets, hardestTwoCMSJets, event_i_));

  // output match inds
  output_match_inds(cms2pfc, "AK52PFC");

  // output jet constituents
  std::vector<bool> pfcJetMask(get_output_mask(pfcJets.size(), cms2pfc));
  output_jets_masked(pfcJets, pfcJetMask, event.pfcs, 'P');

  // handle gen and sim/gen matching
  if (sim_) {
    
    // find gen jets, possibly doing decays
    const VPart * genParticlesPtr(&event.genParticles);
    
    VPart decayedGenParticles;
    if (doDecays_) {

      decayedGenParticles = decay_particles(event.genParticles);
      genParticlesPtr = &decayedGenParticles;

      if (decayedGenParticles.size() == 0) {
        numPythiaFail_++;
        return;
      }
    }

    const VPJ genJets(cluster_gens(*genParticlesPtr));

    // match gen to pfc
    Vint cms2gen(genMatcher_(genJets, hardestTwoCMSJetsCorrected, event_i_));
    std::vector<bool> genJetMask(get_output_mask(genJets.size(), cms2gen));

    // get two final partons
    VPJ allHardPartons(to_pjs(event.hardParticles)), hardPartons;
    hardPartons.push_back(*(allHardPartons.end() - 1));
    hardPartons.push_back(*(allHardPartons.end() - 2));

    // get matched gen to parton
    Vint hard2gen(hardMatcher_(genJets, hardPartons, event_i_));

    // output gen that matched to hard
    for (int i : hard2gen)
      if (i != -1)
        genJetMask[i] = true;

    // gen outputs
    output_match_inds(cms2gen, "AK52GEN");
    output_jets_masked(genJets, genJetMask, *genParticlesPtr, 'G');

    // hard outputs
    output_match_inds(hard2gen, "HARD2GEN");
    for (uint i = 0; i < 2; i++)
      output_particle(hardPartons[i], *(event.hardParticles.end() - 1 - i), 'H', i, PI);
  }

  delete csPFC_;
  delete csGen_;
}

void DijetsAnalyzer::print_progress(std::ostream & ostr) {

  // options
  ostr << "################################################\n"
       << "# OPTIONS\n"
       << "################################################\n"
       << "#\n"
       << "# Pythia          " << PYTHIA_VERSION           << '\n'
       << "# FastJet         " << fastjet::fastjet_version << '\n'
       << "#\n"
       << "# inFilePath      " << inFilePath_  << '\n'
       << "# outFilePath     " << outFilePath_ << '\n'
       << "# fileType        " << fileType_    << '\n'
       << "# num             " << num_         << '\n'
       << "#\n";
  if (sim_) 
  ostr << "# doDecays        " << (doDecays_ ? "true" : "false") << '\n'
       << "#\n";
  ostr << "# jet_R           "  << jet_R_    << '\n'
       << "# minJetPt        "  << minJetPt_ << '\n'
       << "#\n"
       << "# matchCMSPtFrac  " << cmsMatcher_.pt_frac()      << '\n'
       << "# matchCMSDeltaR  " << cmsMatcher_.delta_r_max()  << '\n'
       << "# matchGenPtFrac  " << genMatcher_.pt_frac()      << '\n'
       << "# matchGenDeltaR  " << genMatcher_.delta_r_max()  << '\n'
       << "# matchHardPtFrac " << hardMatcher_.pt_frac()     << '\n'
       << "# matchHardDeltaR " << hardMatcher_.delta_r_max() << '\n';

  // counters
  ostr << '\n'
       << "################################################\n"
       << "# STATS\n"
       << "################################################\n"
       << "# numNoCMSJets              " << numNoCMSJets_   << '\n'
       << "# numPFCsMissingButExpected " << numPFCsMissing_ << '\n';
  if (sim_) 
  ostr << "# numPythiaFail             " << numPythiaFail_ << '\n';

  // jet matcher counts
  ostr << "#\n"
       << "#  " << cmsMatcher_.cols()  << '\n'
       << "#  " << cmsMatcher_.stats() << '\n';
  if (sim_)
  ostr << "#  " << genMatcher_.stats()  << '\n'
       << "#  " << hardMatcher_.stats() << '\n';
  ostr << "#\n"
       << "# " << inFileName_ << ", processed " << event_i_ << " events.\n"
       << '\n';
}

void DijetsAnalyzer::output_event_info(const mod::Event & event) {

  outFile_ << "i "   << event_i_                                       << '\n'
           << "EV "  << event.eventNum                                 << '\n'
           << "NPV " << event.nPV                                      << '\n'
           << "t "   << event.time_s + (double) event.time_us/1000000. << '\n'
           << "RHO " << event.rho                                      << '\n'
           << "iLB " << lb_i_                                          << '\n'
           << "TR ";

  for (uint tr : event.firedTriggers)
    outFile_ << tr << ' ';
  outFile_ << '\n';
}

void DijetsAnalyzer::output_ak5s(const VPJ & pjs, const std::vector<mod::AK5> & ak5s) {
  for (uint i = 0; i < pjs.size(); i++) {
    const fastjet::PseudoJet & pj(pjs[i]);
    const mod::AK5 & ak5(ak5s[pj.user_index()]);

    uint quality(jet_quality(ak5.neutralHadronEnergyFraction,
                             ak5.neutralEmEnergyFraction,
                             ak5.chargedHadronEnergyFraction,
                             ak5.chargedEmEnergyFraction,
                             ak5.nConsts,
                             ak5.nChargedConsts,
                             ak5.eta));

    ak5.write(outFile_, i, pj.phi(), quality);
  }
}

void DijetsAnalyzer::output_jets_masked(const VPJ & jets, const std::vector<bool> & mask, 
                                        const VPart & particles, const char type) {

  assert(jets.size() == mask.size());

  for (uint i = 0; i < jets.size(); i++)
    if (mask[i]) {
      const VPJ jet_consts(jets[i].constituents());
      double jet_phi(jets[i].phi());

      for (auto & pj : jet_consts)
        output_particle(pj, particles[pj.user_index()], type, i, jet_phi);
    }
}

void DijetsAnalyzer::output_particle(const fastjet::PseudoJet & pj, const mod::Particle & p, 
                                            const char type, const int i, const double phi_ref) {

  outFile_ << type                       << ' ' 
           << i                          << ' '
           << pj.pt()                    << ' '
           << pj.rap()                   << ' '
           << phi_fix(pj.phi(), phi_ref) << ' '
           << m_round(pj.m())            << ' '
           << p.pdgid;

  if (p.vertex != -100)
    outFile_ << ' ' << p.vertex;

  outFile_ << '\n';
}

void DijetsAnalyzer::output_match_inds(const Vint & ref2new, const char * names) {

  outFile_ << "Match" << names << ' ';
  for (int i : ref2new) 
    outFile_ << i << ' ';

  outFile_ << '\n';
}

void DijetsAnalyzer::output_lumiblocks() {
  for (uint i = 0; i < lumiBlocks_.size(); i++)
    lumiBlocks_[i].write(outFile_, trimTrigNames_, countTriggers_);
}
