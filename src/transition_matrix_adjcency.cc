#include "transition_matrix.h"
#include "transition_matrix_adjcency.h"

namespace gll {

void TransMatrix::AdjcencyImpl::load(const char* filename)
{
  std::cout << "hi TransMatrix::AdjcencyImpl" << std::endl;
  if (access(filename, F_OK) != 0) {
    std::cerr << "ERROR: can't open graph file." << std::endl;
    exit(1);
  }
  nodesz_ = 0;
  std::ifstream ifs(filename);
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) continue;
    nodesz_++;
    ifs.peek();
  }
  ifs.close();

  mat_ = Matrix(nodesz_ + 1, Array(nodesz_ + 1, 0.0));
  ifs.open(filename);
  unsigned int src_id = 1;
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
