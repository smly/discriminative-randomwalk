#ifndef TRANSITION_MATRIX_ADJCENCY_H
#define TRANSITION_MATRIX_ADJCENCY_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>

namespace gll {

typedef std::vector<double> Array;
typedef std::vector<Array> Matrix;

class TransMatrix::AdjcencyImpl : public TransMatrixImpl {
public:
  void load(const char* filename);
  void normalize();
  double trans(
      unsigned int src,
      unsigned int dst) const;
private:
  Matrix mat_;
  unsigned int nodesz_;
};

}

#endif
