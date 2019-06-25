#ifndef __UMOD_ANALYZER__
#define __UMOD_ANALYZER__

// standard library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace mod {

// typedefs
typedef unsigned int uint;

// structure for representing a lumiblock
struct LumiBlock {

  // constructor from an input stream and a lumiblock number 
  LumiBlock(std::istream & istr, uint lbnum);

  // empty destructor
  ~LumiBlock() {}

  // functions to initialize and write out lumiblock
  void init(std::istream & istr);
  void write(std::ostream & ostr, bool trimNames, bool writeCounts) const;

  // lumiblock info
  uint lumiBlockNum, runNum;

  // sim file will have cross section
  double crossSection;

  // triggers
  std::vector<std::string> triggerNames;
  std::vector<std::pair<uint, uint>> triggerPrescales;
  std::vector<uint> firedCounts;

};

// structure for representing an ak5 jet
struct AK5 {

  // constructor from input stream
  AK5(std::istream & istr);

  // empty destructor
  ~AK5() {}

  // function to write out ak5
  void write(std::ostream & ostr, int i_ = -1, double phi_ = -100, int quality = -100) const;

  // jet index
  uint i;

  // kinematics
  double pt, eta, phi, m;

  // derived properties
  double jec, area;
  uint nConsts, nChargedConsts;
  
  // energy fractions
  double neutralHadronEnergyFraction,
         neutralEmEnergyFraction,
         chargedHadronEnergyFraction,
         chargedEmEnergyFraction;

};

// structure capable of representing a Particle
struct Particle {

  // constructor from input stream
  Particle(std::istream & istr, bool extractVertex = true);

  // constructor from values
  Particle(double px_, double py_, double pz_, double e_, int pdgid_);

  // empty destructor
  ~Particle() {}

  // function to write particle 
  void write(std::ostream & ostr) const;

  // kinematics
  double px, py, pz, e;

  // flavor
  int pdgid;

  // vertexing
  int vertex;

};

// structure for representing an event
struct Event {

  // constructor from input stream, filetype, and event index in file
  Event(std::istream & istr, bool sim, uint i);

  // empty destructor
  ~Event() {}

  // functions to initialize and write out event
  void init(std::istream & istr);
  void write(std::ostream & ostr) const;

  // event type
  const bool sim_;

  // event info
  const uint fileIndexNum;
  uint eventNum, nPV;
  double rho;
  unsigned long long time_s;
  uint time_us;

  // triggers
  std::vector<uint> firedTriggers;

  // AK5s
  std::vector<AK5> ak5s;

  // PFCs
  std::vector<Particle> pfcs;

  // hard particles
  std::vector<Particle> hardParticles;

  // gen particles
  std::vector<Particle> genParticles;

};

// generic analyzer, to be used primarily as a base class
class Analyzer {
public:

  // constructor with input/output filepaths
  Analyzer(const std::string & inFilePath, const std::string & outFilePath, const std::string & fileType, 
           int num, int printEvery, uint outprecision);

  // empty destructor
  ~Analyzer() {}

  // update counts of fired triggers
  void update_fired_trigger_counts(const Event &);

  // parse the file, to be called from the constructor of the derived class
  void parse_file();

  // write the event totals to the output file
  void write_summary();

  // functions to process events and lumiblocks
  virtual void process_event(const Event &) {}
  virtual void process_lumiblock(const LumiBlock &) {}
  virtual void print_progress(std::ostream &) {}

protected:

  // input/output 
  const std::string inFilePath_, outFilePath_, inFileName_, outFileName_, fileType_;
  std::ifstream inFile_;
  std::ofstream outFile_;

  // sim flag
  const bool sim_;

  // precision for the outfile
  const uint outprecision_;

  // number of events to process
  const int num_, printEvery_;

  // lumiblocks
  std::vector<LumiBlock> lumiBlocks_;
  int lb_i_;

  // events
  uint event_i_, numBadEvents_, numGoodEvents_, numTotalEvents_, duration_;

};

} // namespace mod

#endif // __UMOD_ANALYZER__
