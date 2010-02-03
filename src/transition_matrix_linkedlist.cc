#include "transition_matrix.h"
#include "transition_matrix_linkedlist.h"

namespace gll {

void TransMatrix::LinkedListImpl::load(const char* filename)
{
  std::cout << "h8i TransMatrix::LinkedListImpl" << std::endl;
}

void TransMatrix::LinkedListImpl::normalize()
{

}

double TransMatrix::LinkedListImpl::trans(
    unsigned int src, unsigned int dst) const
{
  return 0.0;
}

}
