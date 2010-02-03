#ifndef TRANSITION_MATRIX_LINKEDLIST_H
#define TRANSITION_MATRIX_LINKEDLIST_H

namespace gll {

class TransMatrix::LinkedListImpl : public TransMatrixImpl {
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
