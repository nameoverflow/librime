//
// Copyright RIME Developers
// Distributed under the BSD License
//
// Created by nameoverflow on 2018/11/14.
//

#ifndef RIME_CORRECTOR_H
#define RIME_CORRECTOR_H

#include <rime/common.h>
#include <rime/dict/vocabulary.h>
#include <rime/dict/prism.h>
#include "spelling.h"
#include "algebra.h"

namespace hnswlib {
template<typename dist_t>
class HierarchicalNSW;
template<typename MTYPE>
class SpaceInterface;
};

namespace rime {

class CorrectionCollector {
 public:
  explicit CorrectionCollector(const Syllabary& syllabary): syllabary_(syllabary) {}

  Script Collect(size_t edit_distance);

 private:
  const Syllabary& syllabary_;
};


class Corrector {
 public:
  using Correction = struct {
    float distance;
    SyllableId syllable;
    size_t length;
  };
  RIME_API virtual vector<Correction> ToleranceSearch(const Prism& prism, const string& key, float tolerance) = 0;
};

class EditDistanceCorrector : public Prism,
                              public Corrector {
 public:
  using Distance = uint8_t;

  RIME_API explicit EditDistanceCorrector(const string& file_name) : Prism(file_name) {}

  RIME_API bool Build(const Syllabary& syllabary,
                      const Script* script,
                      uint32_t dict_file_checksum,
                      uint32_t schema_file_checksum) override;

  static Distance LevenshteinDistance(const std::string &s1, const std::string &s2);
  static Distance RestrictedDistance(const std::string& s1, const std::string& s2);

 private:
  vector<Match> SymDeletePrefixSearch(const string& key);
};

class NearSearchCorrector : public Corrector {
 public:
  NearSearchCorrector() = default;

  /// Brute force kNN search
  /// \param prism
  /// \param key
  /// \param tolerance
  /// \return vector of results
  RIME_API vector<Correction> ToleranceSearch(const Prism& prism, const string& key, float tolerance) override;
};

class ANNCorrector : public Corrector,
                     public MappedFile {
 public:
  explicit ANNCorrector(const string& filename): MappedFile(filename) {};
  RIME_API bool Build(const Syllabary& syllabary, Prism &prism);
  RIME_API bool Build(const Script& script, Prism &prism);
  RIME_API bool Load();
  RIME_API bool Save();

  RIME_API vector<Correction> ToleranceSearch(const Prism& prism, const string& key, float tolerance) override;

 private:
  using Metadata = struct {
    size_t data_size = 0;
    size_t dim = 0;
  };
  bool BuildHNSW(Prism &prism, vector<string>& spellings);
  an<hnswlib::HierarchicalNSW<float>> alg_;
  an<hnswlib::SpaceInterface<float>> metric_;
  Metadata *metadata_;
  size_t data_size_ = 0;
  size_t dim_ = 0;
};

} // namespace rime

#endif //RIME_CORRECTOR_H
