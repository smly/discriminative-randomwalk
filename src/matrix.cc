#include "matrix.h"

namespace gll {

void Mat::load(
    const char* graph_fn,
    const char* label_fn,
    const bool symm)
{
  if (access(graph_fn, F_OK) != 0 ||
      access(label_fn, F_OK) != 0) {
    std::cerr << "ERROR: can't open graph or label file." << std::endl;
    exit(1);
  }
  // reading input
  // calc nodesz
  nodesz_ = 0;
  std::ifstream ifs(graph_fn);
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) continue;
    nodesz_++;
    ifs.peek();
  }
  ifs.close();

  // [dummy, 1, .., nodesz]
  mat_ = Matrix(nodesz_+1, Array(nodesz_+1, 0.0)); // **
  ifs.open(graph_fn);
  unsigned int src_id = 1;
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) continue;
    const char* s = line.c_str();
    const unsigned int len = strlen(s);
    unsigned int pos = 0;
    unsigned int elem_num = 0;
    for (unsigned int i = 0; i < len; i++) {
      if (s[i] == ':') elem_num++;
    }
    pos = 0;
    for (unsigned int i = 0; i < elem_num; i++) {
      unsigned int dst_id;
      double w;
      while (pos < len && isspace(s[pos])) pos++;
      dst_id = atoi(s + pos);
      while (pos + 1 < len && s[pos] != ':') pos++;
      w = atof(s + pos + 1);
      while (pos < len && !isspace(s[pos])) pos++;
      mat_[src_id][dst_id] = w;
      //      if (symm) {
      mat_[dst_id][src_id] = w;
      //      }
    }
    ifs.peek();
    src_id++;
  }
  ifs.close();
  ifs.open(label_fn);
  labelsz_ = 0;
  lnodesz_ = 0;
  unodesz_ = 0;
  unsigned int nodesz_ = 0;
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) break;
    const char* s = line.c_str();
    const unsigned int len = strlen(s);
    nodesz_++;
    if (len != 1 || (len >= 1 && isdigit(s[0]))) {
      unsigned int label = atoi(s);
      if (label > labelsz_) {
        labelsz_ = label;
      }
      lnodesz_++;
    } else {
      unodesz_++;
    }
    ifs.peek();
  }
  ifs.close();
  assert(lnodesz_ + unodesz_ == nodesz_);
  lmat_  = LabelMatrix(labelsz_+1, LabelArray());
  lcmat_ = LabelMatrix(labelsz_+1, LabelArray());

  ifs.open(label_fn);
  nodesz_ = 0;
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) break;
    const char* s = line.c_str();
    const unsigned int len = strlen(s);
    nodesz_++;
    if (len != 1 || (len >= 1 && isdigit(s[0]))) {
      unsigned int label = atoi(s);
      lmat_[label].push_back(nodesz_);
    }
  }
  ifs.close();
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    unsigned int indx = 0;
    const unsigned int indx_max = lmat_[lid].size() - 1;
    for (unsigned int nid = 1; nid <= nodesz_; nid++) {
      if (indx <= indx_max && lmat_[lid][indx] == nid) {
        indx++;
      } else {
        lcmat_[lid].push_back(nid);
      }
    }
  }
  // normalized
  normalize();
}

void Mat::normalize() {
  const unsigned int nodesz_ = mat_.size() - 1;
  for (unsigned int i = 1; i <= nodesz_; i++) {
    double sum = 0.0;
    for (unsigned int j = 1; j <= nodesz_; j++) {
      sum += mat_[i][j];
    }
    for (unsigned int j = 1; j <= nodesz_; j++) {
      mat_[i][j] /= sum;
      if (isnan(mat_[i][j])) {
        mat_[i][j] = 0.0;
      }
    }
  }
}

}
