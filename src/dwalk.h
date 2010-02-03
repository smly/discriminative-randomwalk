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

namespace gll {

typedef std::vector<double> Array;
typedef std::vector<Array> Matrix;
typedef std::vector<unsigned int> LabelArray;
typedef std::vector<LabelArray> LabelMatrix;
typedef std::set<unsigned int> NodeSet;

template <class T>
class IncrementalGen {
public:
  IncrementalGen (T t) : curr(t) { }
  T operator () () { return curr++; }
private:
  T curr;
};

struct Alpha {}; // alpha forward variable
struct Beta {}; // beta backward variable
struct Gamma {}; // gamma backward variable

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
  void fixed_algorithm_A(
      const unsigned int bounded_length,
      const std::string& pref_fn,
      const bool map_predict);
  void fixed_algorithm_B(
      const unsigned int bounded_length,
      const std::string& pref_fn,
      const bool map_predict);

  // dwalk_const
  void show_info() const;
  void show_alphabeta(
      const unsigned int t,
      const std::string label,
      const Matrix& mat) const;
  void show_mat(const Matrix& mat) const;
  void show_lmat(const LabelMatrix& lmat) const;
  void show_betweenness(const Matrix& mat) const;
  void show_predict(
      const std::vector<unsigned int>& predict) const;
  const LabelArray get_complement_labeled_nodes (
      const unsigned int label_id) const;
  const LabelArray get_labeled_nodes (
      const unsigned int label_id) const;

private:
  unsigned int nodesz_, lnodesz_, unodesz_;
  unsigned int labelsz_;
  Matrix mat;
  LabelMatrix lmat, lcmat;

  // private func
  void normalize(Matrix& mat);
  void calc_alpha(
      std::vector<Matrix>& alpha_v,
      const unsigned int bounded_length);
  void calc_beta(
      std::vector<Matrix>& beta_v,
      const unsigned int bounded_length);
  void calc_gamma(
      std::vector<Matrix>& beta_v,
      const NodeSet& uset,
      const unsigned int bounded_length);
  void decision(
      const Matrix& b,
      const std::string& pref_fn,
      const bool map_predict_flag);
};

} // end of namespace

#endif
