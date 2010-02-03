#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

namespace gll {

typedef std::vector<double> Array;
typedef std::vector<Array> Matrix;
typedef std::vector<unsigned int> LabelArray;
typedef std::vector<LabelArray> LabelMatrix;

class Mat {
public:
  Mat() { }
  void load(
      const char* graph_fn,
      const char* label_fn,
      const bool symm);
  void normalize();

  unsigned int nodesz_, lnodesz_, unodesz_;
  unsigned int labelsz_;
  Matrix mat_;
  LabelMatrix lmat_, lcmat_;
private:
};

}

#endif
