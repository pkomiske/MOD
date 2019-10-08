#ifndef __DatasetHH__
#define __DatasetHH__

// Standard library headers
#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// HDF5 C++ header
#include "H5Cpp.h"

namespace MOD {

// struct for representing a pfc/gen
struct Particle {

  // allow Particle to be declared without any arguments
  Particle() {}

  // initialize particle from vector iterator
  Particle(std::vector<double>::iterator & it) :
    pt(*it++), y(*it++), phi(*it++), m(*it++),
    pid(int(*it++)), vertex(int(*it++))
  {}

  double pt, y, phi, m;
  int pid, vertex;
};

// base class for representing a mod array containing all the basic operations we want to do to it
class Array {
public:

  Array(H5::H5File h5file, const char * name) :
    dataset_(h5file.openDataSet(name)),
    datatype_(dataset_.getDataType()),
    dataspace_(dataset_.getSpace()),
    numElements_(1)
  {

    // get dimensions
    int ndims_ = dataspace_.getSimpleExtentNdims();
    dims_.resize(ndims_);
    dataspace_.getSimpleExtentDims(dims_.data());

    // store the total number of elements
    for (auto d : dims_) numElements_ *= d;

    // determine chunking
    H5::DSetCreatPropList cproplist(dataset_.getCreatePlist());
    if (H5D_CHUNKED == cproplist.getLayout()) {
      chunks_.resize(ndims_);
      cproplist.getChunk(ndims_, chunks_.data());
    }
  }

  // access methods
  inline H5::DataSet & getDataSet() { return dataset_; }
  inline H5::DataType & getDataType() { return datatype_; }
  inline H5::DataSpace & getDataSpace() { return dataspace_; }
  inline hsize_t getNumElements() const { return numElements_; }
  inline int getNdims() const { return ndims_; }
  inline const std::vector<hsize_t> & getDims() const { return dims_; }
  inline const std::vector<hsize_t> & getChunks() const { return chunks_; }

private:

  // variables used in initialization
  H5::DataSet dataset_;
  H5::DataType datatype_;
  H5::DataSpace dataspace_;
  hsize_t numElements_;

  // variables used in constructor
  int ndims_;
  std::vector<hsize_t> dims_, chunks_;
};

// class representing pfcs_index/gens_index arrays
class IndexArray : public Array {
public:

  IndexArray(H5::H5File h5file, const char * name) :
    Array(h5file, name)
  {
    // read all elements
    dataVec_.resize(getNumElements());
    getDataSet().read(dataVec_.data(), getDataType());
  }

  std::pair<uint32_t, uint32_t> getRangePair(hsize_t i) const {
    return std::make_pair(dataVec_[i], dataVec_[i+1]);
  }

private:
  std::vector<uint32_t> dataVec_;
};

// class representing jets_i/jets_f arrays
template<typename T>
class JetArray : public Array {
public:

  JetArray(H5::H5File h5file, const char * name) :
    Array(h5file, name)
  {
    // get number of columns
    numCols_ = getDims()[1];

    // read all elements
    dataVec_.resize(getNumElements());
    getDataSet().read(dataVec_.data(), getDataType());
  }

  // return row from 2D dataset
  std::vector<T> get2DRow(hsize_t row) const {
    auto start(dataVec_.begin() + row * numCols_);
    return std::vector<T>(start, start + numCols_);
  }

private:
  hsize_t numCols_;
  std::vector<T> dataVec_;
};

// class representing pfcs/gens arrays
class ParticleArray : public Array {
public:

  ParticleArray(H5::H5File h5file, const char * name) :
    Array(h5file, name)
  {
    count_[1] = getDims()[1];
    start_[1] = 0;
  }

  // select hyperslab from dataset dataspace
  std::vector<Particle> get2DRange(const std::pair<uint32_t, uint32_t> & range) {

    // set dataspace to the elements specified by range
    count_[0] = range.second - range.first;
    start_[0] = range.first;
    getDataSpace().selectHyperslab(H5S_SELECT_SET, count_, start_);

    // setup mem dataspace
    H5::DataSpace memDataSpace(2, count_);

    // ensure buffer has enough space
    dataVec_.resize(count_[0] * count_[1]);

    // read
    getDataSet().read(dataVec_.data(), getDataType(), memDataSpace, getDataSpace());

    // make into vector of particles
    std::vector<Particle> particles;
    for (std::vector<double>::iterator it = dataVec_.begin(), end = dataVec_.end(); it != end;)
      particles.emplace_back(it);

    return particles;
  }

private:
  hsize_t count_[2], start_[2];
  std::vector<double> dataVec_;
};

// Dataset class, losely based on Python class of the same name in EnergyFlow package
class Dataset {
public:

  // constructor taking C string as filepath
  Dataset(const char * filepath, const char * collection = "CMS2011AJets") : 
    filepath_(filepath),
    h5file_(filepath_, H5F_ACC_RDONLY),
    collection_(std::string(collection))
  { 
    // make sure we init assuming the correct collection 
    if (collection_ == "CMS2011AJets") initCMS2011AJets(); 
    else throw std::invalid_argument("Collection " + collection_ + " not recognized.");
  }

  // destructor
  ~Dataset() {

    // free arrays
    for (Array* a : arrays_) delete a;
  }

  // common initializations
  void initCMS2011AJets() {

    // set constants
    calledNext_ = false;
    iEvent_ = 0;
    jets_i_loc_ = 0; 
    jets_f_loc_ = 1;
    pfcs_loc_ = gens_loc_ = -1;

    // verify size of double
    static_assert(sizeof(double) == 8, "double must be 8 bytes");

    // setup jets arrays
    arrays_.push_back(new JetArray<int64_t>(h5file_, "jets_i"));
    arrays_.push_back(new JetArray<double>(h5file_, "jets_f"));

    // set num events
    numEvents_ = arrays_[0]->getDims()[0];

    // detect/open pfcs in file
    int next_loc(arrays_.size());
    if (h5file_.nameExists("pfcs")) {
      pfcs_loc_ = next_loc++;
      arrays_.push_back(new ParticleArray(h5file_, "pfcs"));

      pfcs_index_loc_ = next_loc++;
      arrays_.push_back(new IndexArray(h5file_, "pfcs_index"));
    }

    // detect/open gens in file
    if (h5file_.nameExists("gens")) {
      gens_loc_ = next_loc++;
      arrays_.push_back(new ParticleArray(h5file_, "gens"));

      gens_index_loc_ = next_loc++;
      arrays_.push_back(new IndexArray(h5file_, "gens_index"));
    }
  }

  // get index of event in file
  inline hsize_t getIndex() const { return iEvent_; }

  // get total number of events in the file
  inline hsize_t getNumEvents() const { return numEvents_; }

  // determine if we have pfcs
  inline bool hasPFCs() const { return pfcs_loc_ != -1; }

  // determine if we have gens
  inline bool hasGENs() const { return gens_loc_ != -1; }

  // call this to advance each event
  bool next() {

    // do not advance on the very first call
    if (calledNext_) iEvent_++;
    else calledNext_ = true;

    // check if we reached the end
    return getIndex() < getNumEvents();
  }

  // pointer to jets_i array
  JetArray<int64_t>* jets_i_ptr() const {
    return static_cast<JetArray<int64_t>*>(arrays_[jets_i_loc_]);
  }

  // pointer to jets_f array
  JetArray<double>* jets_f_ptr() const {
    return static_cast<JetArray<double>*>(arrays_[jets_f_loc_]);
  }

  // get jets_i corresponding to current index
  std::vector<int64_t> jets_i() const {
    return jets_i_ptr()->get2DRow(getIndex());
  }

  // get jets_f corresponding to current index
  std::vector<double> jets_f() const {
    return jets_f_ptr()->get2DRow(getIndex());
  }

  // generic particle array function
  std::vector<Particle> getParticles(int particles_loc, int index_loc) const {

    // if missing array, return empty vector of particles
    if (particles_loc == -1)
      return std::vector<Particle>();

    // get pfcs range we need
    std::pair<uint32_t, uint32_t> range(static_cast<IndexArray*>(arrays_[index_loc])->getRangePair(getIndex()));

    // return vector of particles
    return static_cast<ParticleArray*>(arrays_[particles_loc])->get2DRange(range);
  }

  // pointer to pfcs array
  ParticleArray* pfcs_ptr() const {
    return static_cast<ParticleArray*>(arrays_[pfcs_loc_]);
  }

  // pointer to gens array
  ParticleArray* gens_ptr() const {
    return static_cast<ParticleArray*>(arrays_[gens_loc_]);
  }

  // get pfcs corresponding to current index
  std::vector<Particle> pfcs() const {
    return getParticles(pfcs_loc_, pfcs_index_loc_);
  }

  // get pfcs corresponding to current index
  std::vector<Particle> gens() const {
    return getParticles(gens_loc_, gens_index_loc_);
  }

private:

  // variables initialized by constructor
  const char * filepath_;
  H5::H5File h5file_;
  const std::string collection_;

  // variables initialized by init
  bool calledNext_;
  hsize_t iEvent_, numEvents_;
  int jets_i_loc_, jets_f_loc_, pfcs_loc_, gens_loc_, pfcs_index_loc_, gens_index_loc_;
  std::vector<Array*> arrays_;
};

} // namespace MOD

#endif // __DatasetHH__