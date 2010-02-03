#include "transition_matrix.h"
#include "transition_matrix_adjcency.h"

namespace gll {

void TransMatrix::AdjcencyImpl::load(const char* filename)
{
  std::cout << "hi TransMatrix::AdjcencyImpl" << std::endl;
}

void TransMatrix::AdjcencyImpl::normalize()
{

}

double TransMatrix::AdjcencyImpl::trans(
    unsigned int src, unsigned int dst) const
{
  return 0.0;
}

}
