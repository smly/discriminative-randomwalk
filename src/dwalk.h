#ifndef DWALK_H
#define DWALK_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "matrix.h"

namespace gll {

typedef std::set<unsigned int> NodeSet;

template <class T>
class IncrementalGen {
public:
  IncrementalGen (T t) : curr(t) { }
  T operator () () { return curr++; }
private:
  T curr;
};

class Dwalk
{
public:
  Dwalk() { }

  void load(
      const char* graph_fn,
      const char* label_fn,
      const bool symm);
  void original(
      const unsigned int bounded_length,
      const std::string& pref_fn,
      const bool denom_flag,
      const bool map_predict);
  void fixedAlgorithmA(
      const unsigned int bounded_length,
      const std::string& pref_fn,
      const bool map_predict);
  void fixedAlgorithmB(
      const unsigned int bounded_length,
      const std::string& pref_fn,
      const bool map_predict);
  unsigned int crossValidation(
      const unsigned int num_fold,
      const unsigned int num_from,
      const unsigned int num_to);

  // dwalk_const
  void showInfo() const;
  void showAlphabeta(
      const unsigned int t,
      const std::string label) const;
  void showMat() const;
  void showLmat() const;
  void showBetweenness() const;
  void showPredict(
      const std::vector<unsigned int>& predict) const;
  const LabelArray getComplementLabeledNodes (
      const unsigned int label_id) const;
  const LabelArray getLabeledNodes (
      const unsigned int label_id) const;

private:
  unsigned int nodesz_, lnodesz_, unodesz_;
  unsigned int labelsz_;
  Matrices m_;
  // alpha forward variables
  std::vector<Matrix> alpha_;
  // beta backward variables
  std::vector<Matrix> beta_;
  // gamma backward variables (not appear in the paper)
  std::vector<Matrix> gamma_;

  // private func
  void calcAlpha(
      const unsigned int bounded_length);
  void calcBeta(
      const unsigned int bounded_length);
  void calcGamma(
      const NodeSet& uset,
      const unsigned int bounded_length);
  void decision(
      const Matrix& b,
      const std::string& pref_fn,
      const bool map_predict_flag);
};

} // end of namespace

#endif
