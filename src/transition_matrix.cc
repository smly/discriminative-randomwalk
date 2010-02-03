#include "transition_matrix.h"
#include "transition_matrix_adjcency.h"
#include "transition_matrix_linkedlist.h"

namespace gll {

TransMatrix::TransMatrix(const transImplType impl_type)
    : type_(impl_type)
{
  switch (type_) {
    case adjcency:
      impl_ = new AdjcencyImpl();
      break;
    case linkedlist:
      impl_ = new LinkedListImpl();
      break;
    default:
      impl_ = NULL;
  }
}

void TransMatrix::load(const char* filename)
{
  if (impl_ != NULL) {
    impl_->load(filename);
  } else {
    // todo: exception
  }
}

} // end of namespace
