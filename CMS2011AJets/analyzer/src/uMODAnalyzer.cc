#include <iomanip>
#include <iostream>

#include "uMODAnalyzer.hh"

namespace mod {


///////////////////////////////////////////////////////////////////////////////
// ANALYZER
///////////////////////////////////////////////////////////////////////////////

Analyzer::Analyzer(const std::string & inFilePath, const std::string & outFilePath, const std::string & fileType, 
                   int num, int printEvery, uint outprecision) :

// filepaths, filenames, files, and options
inFilePath_(inFilePath), outFilePath_(outFilePath), 
inFileName_(inFilePath_.substr(inFilePath_.find_last_of("/") + 1, 36)),
outFileName_(outFilePath_.substr(outFilePath_.find_last_of("/") + 1)),
fileType_(fileType),
inFile_(inFilePath_), outFile_(outFilePath_),
sim_(fileType_ == "sim"),
outprecision_(outprecision),

// counters
num_(num), printEvery_(printEvery), lb_i_(-1), event_i_(0), 
numBadEvents_(0), numGoodEvents_(0), numTotalEvents_(0), duration_(0)
{
  outFile_ << std::setprecision(outprecision_);
}

void Analyzer::parse_file() {

  std::string s, id;
  while (std::getline(inFile_, s)) {

    if (s.empty()) continue;

    std::istringstream line(std::move(s));
    uint val;
    line >> id >> val;

    if (id == "LB") {
      lumiBlocks_.emplace_back(inFile_, val);
      lb_i_++;
      process_lumiblock(lumiBlocks_.back());
    }

    else if (id == "i") {
      Event event(inFile_, sim_, val);

      update_fired_trigger_counts(event);
      process_event(event);
      outFile_ << '\n';
      event_i_++;

      if ((int) event_i_ == num_)
        break;

      if (printEvery_ != -1 && ((int) event_i_ % printEvery_) == 0)
        print_progress(std::cout);
    }

    else if (id == "NumBadEvents")
      numBadEvents_ = val;

    else if (id == "NumGoodEvents")
      numGoodEvents_ = val;

    else if (id == "NumTotalEvents")
      numTotalEvents_ = val;

    else if (id == "Duration(s)")
      duration_ = val;

    else 
      std::cerr << "Unknown line: '" << id << ' ' << val << "'\n";
  }
}

void Analyzer::update_fired_trigger_counts(const Event & event) {

  std::vector<uint> & counts(lumiBlocks_[lb_i_].firedCounts);

  for (uint tr : event.firedTriggers)
    counts[tr]++;
}

void Analyzer::write_summary() {
  outFile_ << "NumBadEvents " << numBadEvents_ << '\n'
           << "NumGoodEvents " << numGoodEvents_ << '\n'
           << "NumTotalEvents " << numTotalEvents_ << '\n'
           << '\n'
           << "Duration(s) " << duration_ << '\n'
           << '\n';
}


///////////////////////////////////////////////////////////////////////////////
// LUMIBLOCK
///////////////////////////////////////////////////////////////////////////////

LumiBlock::LumiBlock(std::istream & istr, uint lbnum) : 
lumiBlockNum(lbnum), crossSection(-1) 
{
  init(istr);
}

void LumiBlock::init(std::istream & istr) {

  std::string s, id;
  while (std::getline(istr, s)) {

    if (s.empty()) {
      firedCounts.resize(triggerNames.size());
      std::fill(firedCounts.begin(), firedCounts.end(), 0);
      return;
    }

    std::istringstream line(std::move(s));
    line >> id;

    if (id == "RN") {
      line >> runNum;
    }

    else if (id == "XS") {
      line >> crossSection;
    }

    else if (id.substr(0, 3) == "HLT") {
      triggerNames.push_back(id);

      uint p1, p2;
      line >> p1 >> p2;
      triggerPrescales.push_back(std::make_pair(p1, p2));
    }

    else 
      std::cerr << "Unknown line in LumiBlock::init starting with '" << id << "'\n";
  }
}

void LumiBlock::write(std::ostream & ostr, bool trimNames, bool writeCounts) const {

  ostr << "LB " << lumiBlockNum << '\n'
       << "RN " << runNum       << '\n';

  if (crossSection >= 0) 
    ostr << "XS " << crossSection << '\n';

  for (uint i = 0, s = triggerNames.size(); i < s; i++) {
    const std::pair<uint, uint> & prescales(triggerPrescales[i]);

    ostr << (trimNames ? triggerNames[i].substr(0, triggerNames[i].find_last_of('_')) 
                       : triggerNames[i])
         << ' ' << prescales.first << ' ' << prescales.second << '\n';
  }

  if (writeCounts) {
    ostr << "FCs ";
    for (uint c : firedCounts)
      ostr << c << ' ';
    ostr << '\n';
  }

  ostr << '\n';
}


///////////////////////////////////////////////////////////////////////////////
// EVENT
///////////////////////////////////////////////////////////////////////////////

Event::Event(std::istream & istr, bool sim, uint i) :
sim_(sim), fileIndexNum(i)
{
  init(istr);
}

void Event::init(std::istream & istr) {

  std::string s, id;
  while (std::getline(istr, s)) {

    if (s.empty())
      return;

    std::istringstream line(std::move(s));
    line >> id;

    if (id == "P")
      pfcs.emplace_back(line);

    else if (sim_ && id == "G")
      genParticles.emplace_back(line, false);

    else if (id == "AK5")
      ak5s.emplace_back(line);

    else if (sim_ && id == "H")
      hardParticles.emplace_back(line, false);

    else if (id == "EV")
      line >> eventNum;

    else if (id == "NPV")
      line >> nPV;

    else if (id == "RHO")
      line >> rho;

    else if (id == "t(s)")
      line >> time_s;

    else if (id == "t(us)")
      line >> time_us;

    else if (id == "TR") {
      uint trig_i;
      while (line >> trig_i && line.good())
        firedTriggers.push_back(trig_i);
    }

    else
      std::cerr << "Unknown line in Event::init starting with '" << id << "'\n";
  }
}

void Event::write(std::ostream & ostr) const {

  ostr << "i "     << fileIndexNum << '\n'
       << "EV "    << eventNum     << '\n'
       << "NPV "   << nPV          << '\n'
       << "t(s) "  << time_s       << '\n'
       << "t(us) " << time_us      << '\n'
       << "RHO "   << rho          << '\n'
       << "TR ";

  for (uint tr : firedTriggers)
    ostr << tr << ' ';
  ostr << '\n';

  for (const AK5 & ak5 : ak5s)
    ak5.write(ostr);

  for (const Particle & pfc : pfcs) {
    ostr << "P ";
    pfc.write(ostr);
  }

  if (sim_) {
    for (const Particle & hard : hardParticles) {
      ostr << "H ";
      hard.write(ostr);
    }

    for (const Particle & gen : genParticles) {
      ostr << "G ";
      gen.write(ostr);
    }
  }

  ostr << '\n';
}


///////////////////////////////////////////////////////////////////////////////
// AK5
///////////////////////////////////////////////////////////////////////////////

AK5::AK5(std::istream & line) {

  line >> i
       >> pt >> eta >> phi >> m
       >> jec >> area
       >> nConsts >> nChargedConsts
       >> neutralHadronEnergyFraction
       >> neutralEmEnergyFraction
       >> chargedHadronEnergyFraction
       >> chargedEmEnergyFraction;
}

void AK5::write(std::ostream & ostr, int i_, double phi_, int quality) const {

  ostr << "AK5 "
       << (i_ == -1 ? i : i_)         << ' '
       << pt                          << ' '
       << eta                         << ' '
       << (phi_ == -100 ? phi : phi_) << ' '
       << m                           << ' '
       << jec                         << ' '
       << area                        << ' '
       << nConsts                     << ' '
       << nChargedConsts              << ' '
       << neutralHadronEnergyFraction << ' '
       << neutralEmEnergyFraction     << ' '
       << chargedHadronEnergyFraction << ' '
       << chargedEmEnergyFraction;

  if (quality != -100)
    ostr << ' ' << quality;

  ostr << '\n';
}


//////////////////////////////////////////////////////////////////////////////
// PARTICLE
///////////////////////////////////////////////////////////////////////////////

Particle::Particle(std::istream & line, bool extractVertex) {

  line >> px >> py >> pz >> e >> pdgid;

  if (extractVertex) 
    line >> vertex;
  else 
    vertex = -2;
}

Particle::Particle(double px_, double py_, double pz_, double e_, int pdgid_) :
px(px_), py(py_), pz(pz_), e(e_), pdgid(pdgid_), vertex(-100)
{}

void Particle::write(std::ostream & ostr) const {

  ostr << px    << ' '
       << py    << ' '
       << pz    << ' '
       << e     << ' '
       << pdgid;

  if (vertex != -100) 
    ostr << ' ' << vertex;

  ostr << '\n';
}

} // end namespace mod
