// C++ std library
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// boost
#include <boost/unordered_map.hpp>

// CMSSW CondFormats
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

// CMSSW DataFormats
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// CMSSW FWCore
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

// CMSSW HLTrigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// CMSSW SimDataFormats
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

// helpful webpages
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMParametersForModules?rev=64
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHighLevelTrigger
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips

#define OUTPRECISION 10

// our subclass of EDProducer
class MODProducer : public edm::EDProducer {
public:

  // typedefs
  typedef unsigned int uint;

  // constructor and destructor
  MODProducer(const edm::ParameterSet &);
  ~MODProducer() {}

  // overloaded functions from EDProducer
  virtual void beginRun(edm::Run &, edm::EventSetup const &);
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endJob();

  // custom functions
  bool validLumiBlock();
  template <class T> bool validateHandle (T & handle, std::string name);

private:

  // variables initialized through Parameter Set
  std::string dataType_, JECPath_, MODFilepath_, logFilepath_;
  std::vector<std::string> certVStr_;
  const uint printEvery_;
  const double minJetPtToWriteJet_, minJetPtToWriteParticles_;

  // constant CMSSW initializations
  edm::InputTag hltConfigTag_, rhoTag_, AK5PFTag_, PFCTag_, primaryVerticesTag_;
  std::string triggerCategory_;
  HLTConfigProvider hltConfig_;

  // files
  std::ofstream logFile_, MODFile_;

  // stringstream
  std::ostringstream MODoss_;

  // data structure for validating good runs
  boost::unordered_map<uint, std::vector<std::pair<uint, uint> > > certMap_;

  // if we have generated particles
  const bool haveData_;

  // integers
  uint iEvent_, iRun_, badEventCount_, runNum_, lumiBlock_, eventNum_;

  // timing
  time_t startTime_;

  // needs to be initilized in constructor
  FactorizedJetCorrector * AK5JetCorrector_;

  // trigger names
  const std::vector<std::string> * triggerNames_;
  edm::ParameterSetID parameterSetID_;
  bool firstEventInRun_;

  // cross section
  double crossSection_;
};

MODProducer::MODProducer(const edm::ParameterSet & pSet) :

// options contained in the ParameterSet
dataType_(pSet.getParameter<std::string>("data_type")),
JECPath_(pSet.getParameter<std::string>("jec_path")),
MODFilepath_(pSet.getParameter<std::string>("mod_filepath")),
logFilepath_(pSet.getParameter<std::string>("log_filepath")),
certVStr_(pSet.getParameter<std::vector<std::string> >("cert_str")),
printEvery_(pSet.getParameter<uint>("print_every")),
minJetPtToWriteJet_(pSet.getParameter<double>("min_jetpt_jet")),
minJetPtToWriteParticles_(pSet.getParameter<double>("min_jetpt_particles")),

// CMSSW constants
hltConfigTag_("TriggerResults", "", "HLT"),
rhoTag_("kt6PFJets", "rho"),
AK5PFTag_("ak5PFJets"),
PFCTag_("particleFlow"),
primaryVerticesTag_("offlinePrimaryVertices"),
triggerCategory_("Jet"),
hltConfig_(),

// our initializations
logFile_(logFilepath_.c_str(), std::ios_base::out | std::ios_base::app),
MODFile_(MODFilepath_.c_str()),
MODoss_(""),
certMap_(),
haveData_(dataType_ == "cms"),
iEvent_(0),
iRun_(0),
badEventCount_(0),
startTime_(time(NULL))
{
  // parse certVStr into memory if using data
  // note that this scheme carefully coordinates with the python script
  if (haveData_) {
    for (uint i = 0; i < certVStr_.size(); i++) {
      std::string & line(certVStr_[i]);
      
      // get location of colon
      size_t colon_pos(line.find(':'));
      size_t begin(colon_pos + 1), end;

      // extract run number from before the colon in the line
      uint run_number(std::stoi(line.substr(0, colon_pos)));

      // get valid lumi ranges
      std::vector<std::pair<uint, uint> > lumi_ranges;
      while ((end = line.find(',', begin)) != std::string::npos) {
        std::istringstream iss(line.substr(begin, end - begin));
        
        uint lumi_begin, lumi_end;
        iss >> lumi_begin >> lumi_end;

        lumi_ranges.push_back(std::make_pair(lumi_begin, lumi_end));
        begin = end + 1;
      }

      // insert key, value pair into certMap_
      certMap_[run_number] = lumi_ranges;
    }
  }

  // initialize jet corrector
  std::vector<JetCorrectorParameters> JCPvec;
  JCPvec.push_back(JetCorrectorParameters(JECPath_ + "L1FastJet_AK5PF.txt"));
  JCPvec.push_back(JetCorrectorParameters(JECPath_ + "L2Relative_AK5PF.txt"));
  JCPvec.push_back(JetCorrectorParameters(JECPath_ + "L3Absolute_AK5PF.txt"));
  JCPvec.push_back(JetCorrectorParameters(JECPath_ + "L2L3Residual_AK5PF.txt"));
  AK5JetCorrector_ = new FactorizedJetCorrector(JCPvec);

  // output some information to the logFile
  logFile_ << "===========================================================\n"
           << "----------------------- MODProducer -----------------------\n"
           << "===========================================================\n\n"
           << "Data Type:    " << dataType_ << '\n'
           << "MOD Filepath: " << MODFilepath_ << '\n'
           << "JEC Path:     " << JECPath_ << '\n'
           << std::flush;
  logFile_ << std::setprecision(OUTPRECISION);

  // outprecision for MOD file
  MODFile_ << std::setprecision(OUTPRECISION);
  MODoss_  << std::setprecision(OUTPRECISION);   
}

// executes at beginning of run
void MODProducer::beginRun(edm::Run & run, edm::EventSetup const & eventSetup) {
  bool changed(true);
  std::string process(hltConfigTag_.process());

  // initialize trigger configuration (may change per run)
  if (!hltConfig_.init(run, eventSetup, process, changed)) {
    std::stringstream message("beginRun ");
    message << iRun_ << ": ERROR - hltConfig failed to initialize with process" << process << '\n';
    std::cerr << message.str() << MODFilepath_ << '\n';
    logFile_ << message.str();
  }
  else logFile_ << "beginRun " << iRun_ 
                << ": hltConfig successful init, process " << process 
                << ", changed " << changed << '\n';

  if (!haveData_) {
    edm::Handle<GenRunInfoProduct> genRunInfo;
    run.getByLabel("generator", genRunInfo);
    crossSection_ = genRunInfo->crossSection();
  }
  
  // get trigger names
  triggerNames_ = & hltConfig_.datasetContent(triggerCategory_);
  firstEventInRun_ = true;

  // update run number
  iRun_++;
}

// checks if internal runNum and lumiBlock are valid according to certMap
bool MODProducer::validLumiBlock() {
  std::vector<std::pair<uint, uint> > & validRanges(certMap_[runNum_]);
  for (uint i = 0; i < validRanges.size(); i++) {
    std::pair<uint, uint> & range(validRanges[i]);
    if (range.first <= lumiBlock_ && lumiBlock_ <= range.second) return true;
  }
  return false;
}

// ensure that handle is valid and write error messages if not
template <class T> bool MODProducer::validateHandle(T & handle, std::string name) {
  if (!handle.isValid()) {
    std::ostringstream message("ERROR - invalid ");
    message << name << ", (MODFilepath runNum lumiBlock eventNum) " << MODFilepath_ << ' ' 
                                                                    << runNum_      << ' ' 
                                                                    << lumiBlock_   << ' ' 
                                                                    << eventNum_    << '\n';
    std::cerr <<  message.str();
    logFile_ << message.str();
    return false;
  }
  return true;
}

// get primary vertex of candidate (no need for this to be in the class since it doesn't use any data members)
int getPrimaryVertex(const reco::VertexCollection & vertices, const reco::PFCandidate & pfc, bool useTrack = true) {
  int iVertex(-1);
  if (useTrack) {
    reco::TrackRef const & track(pfc.trackRef());

    // check for null track, indicating neutral particle (given vertex 0)
    if (track.isNull()) return 0;

    float maxweight(0);
    for (uint i = 0; i < vertices.size(); i++) {
      float weight = vertices[i].trackWeight(track);
      if (weight > maxweight) {
        maxweight = weight;
        iVertex = i;
      }
    }
    if (iVertex != -1) return iVertex;
  }

  // find closest vertex in z
  double dzmin(1000000), pfc_z(pfc.vz());
  bool foundVertex(false);
  for (uint i = 0; i < vertices.size(); i++) {
    double dz = fabs(pfc_z - vertices[i].z());
    if (dz < dzmin) {
      dzmin = dz;
      iVertex = i;
      foundVertex = true;
    }
  }
  if (foundVertex) return iVertex;

  // return -1 as default behavior if everything fails
  return -1;
}

void MODProducer::produce(edm::Event & event, const edm::EventSetup & eventSetup) {

  // get event properties
  runNum_ = event.id().run();
  eventNum_ = event.id().event();
  uint lumiBlock(event.luminosityBlock());

  // determine if we should write full trigger info
  bool writeFullTriggers(false);
  if (firstEventInRun_ || lumiBlock != lumiBlock_) writeFullTriggers = true;
  lumiBlock_ = lumiBlock;

  // check that event is valid, skip if not
  if (haveData_ && !validLumiBlock()) {
    badEventCount_++;
    return;
  }

  // get primary vertices
  edm::Handle<reco::VertexCollection> primaryVertices;
  event.getByLabel(primaryVerticesTag_, primaryVertices);
  if (!validateHandle(primaryVertices, "primary vertex collection")) return;
  
  // get trigger results
  edm::TriggerResultsByName triggers(event.triggerResultsByName(hltConfigTag_.process()));
  edm::ParameterSetID pSetID = triggers.parameterSetID();

  // store parameterSetID for comparison with later runs
  if (firstEventInRun_) {
    parameterSetID_ = pSetID;
    firstEventInRun_ = false;
  }

  // check that triggers have not change unexpectedly
  else if (pSetID != parameterSetID_) {
    logFile_ << "ERROR: parameterSetID does not match expectation, triggers appear to be changing mid-run! "
             << "(RN, EV, LB) " << runNum_ << ' ' << eventNum_ << ' ' << lumiBlock_ << '\n';
    std::cerr << "ERROR! Unexpected trigger behavior, check log file " << MODFilepath_ << '\n';
  }

  // output all triggers if this is a new lumiBlock
  if (writeFullTriggers) {
    MODoss_ << "LB " << lumiBlock_ << '\n'
             << "RN " << runNum_    << '\n';

    // output cross section if we have gen
    if (!haveData_)
      MODoss_ << "XS " << crossSection_ << '\n';

    // iterate over triggers and output prescales
    for (uint i = 0; i < triggerNames_->size(); i++) {
      const std::string & trigName((*triggerNames_)[i]);
      const std::pair<int, int> & prescales(hltConfig_.prescaleValues(event, eventSetup, trigName));
      MODoss_ << trigName         << ' '
               << prescales.first  << ' '
               << prescales.second << '\n';
    }
    MODoss_ << '\n';
  }

  // get rho estimation
  edm::Handle<double> rhoHandle;
  event.getByLabel(rhoTag_, rhoHandle);
  double rho(*rhoHandle);

  // event wide output
  MODoss_ << "i "      << iEvent_                          << '\n'
           << "EV "    << eventNum_                        << '\n'
           << "NPV "   << primaryVertices->size()          << '\n'
           << "t(s) "  << event.time().unixTime()          << '\n'
           << "t(us) " << event.time().microsecondOffset() << '\n'
           << "RHO "   << rho                              << '\n';

  // iterate over triggers to output which ones fired
  MODoss_ << "TR ";
  for (uint i = 0; i < triggerNames_->size(); i++)
    if (triggers.accept((*triggerNames_)[i])) MODoss_ << i << ' ';
  MODoss_ << '\n';

  // get AK5 jets
  edm::Handle<reco::PFJetCollection> AK5s;
  event.getByLabel(AK5PFTag_, AK5s);
  if (!validateHandle(AK5s, "AK5 collection")) return;

  // iterate over jets
  bool writeParticles(false);
  uint ak5_i(0);
  for (reco::PFJetCollection::const_iterator it = AK5s->begin(), end = AK5s->end(); it != end; it++, ak5_i++) {
    double pt(it->pt()), eta(it->eta()), area(it->jetArea());

    // get JEC factor
    AK5JetCorrector_->setJetPt(pt);
    AK5JetCorrector_->setJetEta(eta);
    AK5JetCorrector_->setJetA(area);
    AK5JetCorrector_->setRho(rho);
    double jec(AK5JetCorrector_->getCorrection());
    double pt_with_jec(pt*jec);

    // determine if we should write the particles later
    if (pt_with_jec >= minJetPtToWriteParticles_) writeParticles = true;

    // determine if jet passes cut to be written out
    if (pt_with_jec >= minJetPtToWriteJet_ || ak5_i < 2) { 
      MODoss_ << "AK5 " 
              << ak5_i                             << ' '
              << pt                                << ' '
              << eta                               << ' '
              << it->phi()                         << ' '
              << it->mass()                        << ' '
              << jec                               << ' '
              << area                              << ' '
              << it->nConstituents()               << ' '
              << it->chargedMultiplicity()         << ' '
              << it->neutralHadronEnergyFraction() << ' '
              << it->neutralEmEnergyFraction()     << ' '
              << it->chargedHadronEnergyFraction() << ' '
              << it->chargedEmEnergyFraction()     << '\n';
    }
  }

  // write down particles in event
  if (writeParticles) {

    // get PFCs from event
    edm::Handle<reco::PFCandidateCollection> PFCs;
    event.getByLabel(PFCTag_, PFCs);
    if (!validateHandle(PFCs, "PFC collection")) return;

    // iterate over PFCs
    for(reco::PFCandidateCollection::const_iterator it = PFCs->begin(), end = PFCs->end(); it != end; it++) {
      int vertex(getPrimaryVertex(*primaryVertices, *it));

      MODoss_ << "P "
              << it->px()     << ' '
              << it->py()     << ' '
              << it->pz()     << ' '
              << it->energy() << ' '
              << it->pdgId()  << ' '
              << vertex       << '\n';
    }

    // write down gen particles
    if (!haveData_) {

      // obtain gen particles from event
      edm::Handle<reco::GenParticleCollection> genParticles;
      event.getByLabel("genParticles", genParticles);

      // iterate over gen candidates
      for(reco::GenParticleCollection::const_iterator it = genParticles->begin(), end = genParticles->end(); it != end; it++) {

        // output gen particle
        if (it->status() == 1)
          MODoss_ << "G "
                  << it->px()     << ' '
                  << it->py()     << ' '
                  << it->pz()     << ' '
                  << it->energy() << ' '
                  << it->pdgId()  << '\n';

        // output hard process particle
        else if (it->status() == 3)
          MODoss_ << "H "
                  << it->px()     << ' '
                  << it->py()     << ' '
                  << it->pz()     << ' '
                  << it->energy() << ' '
                  << it->pdgId()  << '\n';
      }
    }
  }

  // newline between events
  MODoss_ << '\n';

  // increase event counter
  iEvent_++;

  // write to files
  if (iEvent_ % printEvery_ == 0) {
    logFile_ << iEvent_ << " valid events processed in " << difftime(time(NULL), startTime_) << "s\n" << std::flush;
    MODFile_ << MODoss_.str() << std::flush;

    // reset MODoss
    MODoss_.str("");
    MODoss_ << std::setprecision(OUTPRECISION);
  }
}

// final output information
void MODProducer::endJob() {
  double duration(difftime(time(NULL), startTime_));
  std::ostringstream oss;

  oss << "NumBadEvents "   << badEventCount_             << '\n'
      << "NumGoodEvents "  << iEvent_                    << '\n'
      << "NumTotalEvents " << (badEventCount_ + iEvent_) << "\n\n"
      << "Duration(s) "    << duration                   << "\n\n";

  // output remaining MODoss as well as total information
  MODFile_ << MODoss_.str() << oss.str();
  logFile_ << oss.str();
}

// declare plugin to cmssw
DEFINE_FWK_MODULE(MODProducer);
