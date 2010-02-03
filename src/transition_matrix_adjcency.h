#ifndef TRANSITION_MATRIX_ADJCENCY_H
#define TRANSITION_MATRIX_ADJCENCY_H

#include <vector>

namespace gll {

class TransMatrix::AdjcencyImpl : public TransMatrixImpl {
public:
  void load(const char* filename);
  void normalize();
  double trans(
      unsigned int src,
      unsigned int dst) const;
private:
};

}

#endif
