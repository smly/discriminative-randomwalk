#include "dwalk.h"

void
Dwalk::calc_alpha(
    std::vector<Matrix>& alpha_v,
    const Matrix& mat,
    const LabelMatrix& lmat,
    const LabelMatrix& lcmat,
    const unsigned int L,
    const unsigned int labelsz,
    const unsigned int nodesz)
{
  Matrix alpha(labelsz+1, Array(nodesz+1)); // **
  Matrix alpha_next(labelsz+1, Array(nodesz+1)); // **

  // alpha first
  std::cout << "calculate alpha variables: ";
  for (unsigned int lid = 1; lid <= labelsz; lid++) {
    for (unsigned int sid = 1; sid <= nodesz; sid++) {
      const unsigned int nysz = lmat[lid].size();
      double sum = 0.0;
      for (unsigned int i = 0; i < lmat[lid].size(); i++) {
        const unsigned int did = lmat[lid][i];
        sum += (1.0 / nysz) * mat[did][sid];
      }
      //      std::cout << lid << "," << sid << "," << "nysz=" << nysz << "," << sum << std::endl;
      alpha[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  alpha_v[1] = alpha;
  for (unsigned int iter = 2; iter <= L; iter++) {
    // alpha next
    for (unsigned int lid = 1; lid <= labelsz; lid++) {
      for (unsigned int sid = 1; sid <= nodesz; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < lcmat[lid].size(); i++) {
          const unsigned int did = lcmat[lid][i];
          sum += alpha[lid][did] * mat[did][sid];
        }
        alpha_next[lid][sid] = sum;
      }
    }
    alpha = alpha_next;
    alpha_v[iter] = alpha;
    std::cout << "." << std::flush;
  }
  std::cout << "done" << std::endl;
}

void
Dwalk::calc_beta(
    std::vector<Matrix>& beta_v,
    const Matrix& mat,
    const LabelMatrix& lmat,
    const LabelMatrix& lcmat,
    const unsigned int L,
    const unsigned int labelsz,
    const unsigned int nodesz)
{
  // calculate backward variable beta
  Matrix beta(labelsz+1, Array(nodesz+1));
  Matrix beta_next(labelsz+1, Array(nodesz+1));
  // beta first
  std::cout << "calculate beta variables:  " << std::flush;
  for (unsigned lid = 1; lid <= labelsz; lid++) {
    for (unsigned sid = 1; sid <= nodesz; sid++) {
      double sum = 0.0;
      for (unsigned int i = 0; i < lmat[lid].size(); i++) {
        const unsigned int did = lmat[lid][i];
        sum += mat[sid][did];
      }
      beta[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  beta_v[1] = beta;
  // beta next
  for (unsigned int iter = 2; iter <= L; iter++) {
    for (unsigned int lid = 1; lid <= labelsz; lid++) {
      for (unsigned int sid = 1; sid <= nodesz; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < lcmat[lid].size(); i++) {
          const unsigned int did = lcmat[lid][i];
          sum += beta[lid][did] * mat[sid][did];
        }
        beta_next[lid][sid] = sum;
      }
    }
    beta = beta_next;
    beta_v[iter] = beta;
    std::cout << "." << std::flush;
  }
  std::cout << "done" << std::endl;
}

void
Dwalk::normalize(Matrix& mat)
{
  const unsigned int nodesz = mat.size() - 1;
  for (unsigned int i = 1; i <= nodesz; i++) {
    double sum = 0.0;
    for (unsigned int j = 1; j <= nodesz; j++) {
      sum += mat[i][j];
    }
    for (unsigned int j = 1; j <= nodesz; j++) {
      mat[i][j] /= sum;
    }
  }
}

void
Dwalk::show_predict(const std::vector<unsigned int>& predict)
{
  const unsigned int nodesz = predict.size() - 1;
  for (unsigned int i = 1; i <= nodesz; i++) {
    std::cout << predict[i] << ' ';
  }
  std::cout << std::endl;
}

void
Dwalk::show_alphabeta(
    const unsigned int t,
    const std::string label,
    const Matrix& mat)
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(3);
  const unsigned int nodesz = mat[0].size() - 1;
  const unsigned int labelsz = mat.size() - 1;
  std::cout << "forward-backward var:  "
            << label << std::endl
            << "iter:                  "
            << t << std::endl;
  std::cout << "Number of nodes:       "
            << nodesz << std::endl;
  std::cout << "Number of labels:      "
            << labelsz << std::endl;
  for (unsigned int j = 1; j <= labelsz; j++) {
    for (unsigned int i = 1; i <= nodesz; i++) {
      std::cout << mat[j][i] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Dwalk::show_lmat(const LabelMatrix& lmat)
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(3);
  const unsigned int labelsz = lmat.size() - 1;
  // show lmat
  std::cout << "Number of labels:      "
            << labelsz << std::endl;
  for (unsigned int i = 1; i <= labelsz; i++) {
    std::cout << "label(" << i << "): ";
    for (unsigned int j = 0; j < lmat[i].size(); j++) {
      std::cout << lmat[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Dwalk::show_betweenness(const Matrix& mat)
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(3);
  const unsigned int labelsz = mat.size() - 1;
  const unsigned int nodesz = mat[0].size() - 1;
  for (unsigned int i = 1; i <=labelsz; i++) {
    for (unsigned int j = 1; j <= nodesz; j++) {
      std::cout << mat[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Dwalk::show_mat(const Matrix& mat)
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(3);
  // show mat
  for (unsigned int i = 1; i <= nodesz; i++) {
    for (unsigned int j = 1; j <= nodesz; j++) {
      std::cout << mat[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Dwalk::load(
    const char* graph_fn,
    const char* label_fn)
{
  if (access(graph_fn, F_OK) != 0 ||
      access(label_fn, F_OK) != 0) {
    std::cerr << "ERROR: can't open graph or label file." << std::endl;
    exit(1);
  }
  // reading input
  // calc nodesz
  nodesz = 0;
  std::ifstream ifs(graph_fn);
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) continue;
    nodesz++;
    ifs.peek();
  }
  ifs.close();

  // [dummy, 1, .., nodesz]
  mat = Matrix(nodesz+1, Array(nodesz+1, 0.0)); // **
  ifs.open(graph_fn);
  unsigned int src_id = 1;
  while (!ifs.eof()) {
    namespace al = boost::algorithm;
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
      mat[src_id][dst_id] = w;
      mat[dst_id][src_id] = w;
    }
    ifs.peek();
    src_id++;
  }
  ifs.close();
  ifs.open(label_fn);
  labelsz = 0;
  unsigned int nodesz_ = 0;
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) break;
    const char* s = line.c_str();
    const unsigned int len = strlen(s);
    nodesz_++;
    if (len != 1 || (len >= 1 && s[0] != '.')) {
      unsigned int label = atoi(s);
      if (label > labelsz) labelsz = label;
    }
    ifs.peek();
  }
  ifs.close();
  assert(nodesz == nodesz_);
  lmat  = LabelMatrix(labelsz+1, LabelArray());
  lcmat = LabelMatrix(labelsz+1, LabelArray());

  ifs.open(label_fn);
  nodesz_ = 0;
  while (!ifs.eof()) {
    std::string line;
    getline(ifs, line);
    if (line.empty()) break;
    const char* s = line.c_str();
    const unsigned int len = strlen(s);
    nodesz_++;
    if (len != 1 || (len >= 1 && s[0] != '.')) {
      unsigned int label = atoi(s);
      lmat[label].push_back(nodesz_);
    }
  }
  ifs.close();
  for (unsigned int lid = 1; lid <= labelsz; lid++) {
    unsigned int indx = 0;
    const unsigned int indx_max = lmat[lid].size() - 1;
    for (unsigned int nid = 1; nid <= nodesz; nid++) {
      if (indx <= indx_max && lmat[lid][indx] == nid) {
        indx++;
      } else {
        lcmat[lid].push_back(nid);
      }
    }
  }

  // normalized
  normalize(mat);
}

void
Dwalk::show_info()
{
  std::cout << "Number of nodes:           "
            << nodesz << std::endl
            << "Number of labels:          "
            << lmat.size() - 1 << std::endl
            << "Number of labeled nodes:   "
            << labelsz << std::endl
            << "Number of unlabeled nodes: "
            << nodesz - labelsz << std::endl;

  // debug
  /*
  std::cout << "[> transition matrix" << std::endl;
  show_mat(mat);
  std::cout << "[> labeled node" << std::endl;
  show_lmat(lmat);
  std::cout << "[> complement labeled node" << std::endl;
  show_lmat(lcmat);
  */
}
LabelArray
Dwalk::get_complement_labeled_nodes (const unsigned int label_id)
{
  return lcmat[label_id];
}

LabelArray
Dwalk::get_labeled_nodes (const unsigned int label_id)
{
  return lmat[label_id];
}
void
Dwalk::go(
    const unsigned int L,
    const std::string& pref_fn)
{
  // calculate forward variable alpha
  // labelsz, nodesz, mat, lmat
  std::vector<Matrix> alpha_all(L+1);
  std::vector<Matrix> beta_all(L+1);
  calc_alpha(alpha_all, mat, lmat, lcmat, L, labelsz, nodesz);
  //  show_alphabeta(L, "alpha", alpha_all[L]);
  calc_beta(beta_all, mat, lmat, lcmat, L, labelsz, nodesz);
  //  show_alphabeta(L, "beta", beta_all[L]);
  /*
  // show mat
  show_lmat(lmat);
  show_lmat(lcmat);
  show_mat(mat);
  */
  // calc bounded dwalks betweeness
  Matrix b(labelsz+1, Array(nodesz+1));
  for (unsigned i = 1; i <= labelsz; i++) {
    double denom = 0;
    for (unsigned l = 1; l <= L; l++) {
      for (unsigned j = 0; j < lmat[i].size(); j++) {
        const unsigned int node_id = lcmat[i][j];
        denom += alpha_all[l][i][node_id];
      }
    }
    b[i][0] = 1.0 / denom;
    for (unsigned int j = 1; j <= nodesz; j++) {
      b[i][j] = b[i][0];
    }
  }
  for (unsigned int i = 1; i <= labelsz; i++) {
    for (unsigned int j = 1; j <= nodesz; j++) {
      double sum = 0;
      for (unsigned int l = 2; l <= L; l++) {
        for (unsigned int t = 1; t <= l - 1; t++) {
          sum += alpha_all[l][i][j] * beta_all[l-t][i][j];
        }
      }
      b[i][j] *= sum;
    }
  }
  // MAP decision rule
  Matrix map(labelsz+1, Array(nodesz+1));
  Array prob_y(labelsz+1);
  std::vector<unsigned int> predict(nodesz+1);

  unsigned int labeled_nodesz = 0;
  for (unsigned int i = 1; i <= nodesz; i++) {
    double sum = 0.0;
    for (unsigned int j = 1; j <= labelsz; j++) {
      sum += b[j][i];
    }
    map[0][i] = sum;
  }
  for (unsigned int i = 1; i <= nodesz; i++) {
    for (unsigned int j = 1; j <= labelsz; j++) {
      map[j][i] = b[j][i] / map[0][i];
    }
  }
  for (unsigned int i = 1; i <= labelsz; i++) {
    labeled_nodesz += lmat[i].size();
  }
  for (unsigned int i = 1; i <= labelsz; i++) {
    prob_y[i] = static_cast<double>(lmat[i].size()) / labeled_nodesz;
    std::cout << "prob_y[" << i << "] = " << prob_y[i] << std::endl;
  }
  for (unsigned int i = 1; i <= nodesz; i++) {
    unsigned int predict_tmp = -1;
    double argmax = -0.1;
    for (unsigned int j = 1; j <= labelsz; j++) {
      double tmp = map[j][i] * prob_y[j];
      if (argmax < tmp) {
        predict_tmp = j;
        argmax = tmp;
      }
    }
    predict[i] = predict_tmp;
  }

  // todo
  // saving pref_fn + name
  std::cout << "[> betweenness:" << std::endl;
  show_betweenness(b);
  std::cout << "[> normalized:" << std::endl;
  show_betweenness(map);
  std::cout << "[> classfied result using a MAP:" << std::endl;
  show_predict(predict);
}
