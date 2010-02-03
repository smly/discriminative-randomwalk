#ifndef GRAPH_INTERFACE_H
#define GRAPH_INTERFACE_H

#include <iostream>

namespace gll {

enum transImplType {
  adjcency = 0,
  linkedlist,
};

class TransMatrixImpl;

// transition matrix
class TransMatrix {
public:
  TransMatrix(const transImplType impl_type);
  ~TransMatrix() { }

  void load(const char* filename);

private:
  transImplType type_;
  // transition matrix interface
  TransMatrixImpl * impl_;
  // transition matrix implementation
  class AdjcencyImpl;
  class LinkedListImpl;
};

class TransMatrixImpl {
public:
  virtual void load(const char* filename) = 0;
  virtual void normalize() = 0;
  virtual double trans(
      unsigned int src,
      unsigned int dst) const = 0;
private:
};

}

#endif
