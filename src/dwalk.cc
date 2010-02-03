#include "dwalk.h"

namespace gll {

void
Dwalk::calc_alpha(
    std::vector<Matrix>& alpha_v,
    const unsigned int bounded_length)
{
  Matrix alpha(labelsz_ + 1, Array(nodesz_ + 1, 0.0));
  Matrix alpha_next(labelsz_ + 1, Array(nodesz_ + 1, 0.0));
  // calc \alpha^y(q,1)
  std::cout << "calculate alpha variables: ";
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    for (unsigned int sid = 1; sid <= nodesz_; sid++) {
      const unsigned int nysz = lmat[lid].size();
      double sum = 0.0;
      for (unsigned int i = 0; i < lmat[lid].size(); i++) {
        const unsigned int did = lmat[lid][i];
        sum += (1.0 / nysz) * mat[did][sid];
        assert(isnan(sum) == 0);
      }
      alpha[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  alpha_v[1] = alpha;
  for (unsigned int iter = 2; iter <= bounded_length; iter++) {
    for (unsigned int lid = 1; lid <= labelsz_; lid++) {
      for (unsigned int sid = 1; sid <= nodesz_; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < lcmat[lid].size(); i++) {
          const unsigned int did = lcmat[lid][i];
          sum += alpha[lid][did] * mat[did][sid];
          assert(isnan(sum) == 0);
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
    const unsigned int bounded_length)
{
  Matrix beta(labelsz_+1, Array(nodesz_+1, 0.0));
  Matrix beta_next(labelsz_+1, Array(nodesz_+1, 0.0));
  // calc \beta^y(q,1)
  std::cout << "calculate beta variables:  " << std::flush;
  for (unsigned lid = 1; lid <= labelsz_; lid++) {
    for (unsigned sid = 1; sid <= nodesz_; sid++) {
      double sum = 0.0;
      for (unsigned int i = 0; i < lmat[lid].size(); i++) {
        const unsigned int did = lmat[lid][i];
        sum += mat[sid][did];
        assert(isnan(sum) == 0);
      }
      beta[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  beta_v[1] = beta;
  for (unsigned int iter = 2; iter <= bounded_length; iter++) {
    for (unsigned int lid = 1; lid <= labelsz_; lid++) {
      for (unsigned int sid = 1; sid <= nodesz_; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < lcmat[lid].size(); i++) {
          const unsigned int did = lcmat[lid][i];
          sum += beta[lid][did] * mat[sid][did];
          assert(isnan(sum) == 0);
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
Dwalk::calc_gamma(
    std::vector<Matrix>& gamma_v,
    const NodeSet& uset,
    const unsigned int bounded_length)
{
  Matrix gamma(labelsz_+1, Array(nodesz_+1, 0.0));
  Matrix gamma_next(labelsz_+1, Array(nodesz_+1, 0.0));
  // calc \gamma^y(q,1)
  std::cout << "calculate gamma variables: " << std::flush;
  for (unsigned lid = 1; lid <= labelsz_; lid++) {
    for (unsigned sid = 1; sid <= nodesz_; sid++) {
      double sum = 0.0;
      for (NodeSet::iterator it = uset.begin(); it != uset.end(); it++) {
        const unsigned int did = *it;
        sum += mat[sid][did];
        assert(isnan(sum) == 0);
      }
      gamma[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  gamma_v[1] = gamma;
  for (unsigned int iter = 2; iter <= bounded_length; iter++) {
    for (unsigned int lid = 1; lid <= labelsz_; lid++) {
      for (unsigned int sid = 1; sid <= nodesz_; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < lcmat[lid].size(); i++) {
          const unsigned int did = lcmat[lid][i];
          sum += gamma[lid][did] * mat[sid][did];
          assert(isnan(sum) == 0);
        }
        gamma_next[lid][sid] = sum;
      }
    }
    gamma = gamma_next;
    gamma_v[iter] = gamma;
    std::cout << "." << std::flush;
  }
  std::cout << "done" << std::endl;
}

void
Dwalk::normalize(Matrix& mat)
{
  const unsigned int nodesz_ = mat.size() - 1;
  for (unsigned int i = 1; i <= nodesz_; i++) {
    double sum = 0.0;
    for (unsigned int j = 1; j <= nodesz_; j++) {
      sum += mat[i][j];
    }
    for (unsigned int j = 1; j <= nodesz_; j++) {
      mat[i][j] /= sum;
      if (isnan(mat[i][j])) {
        mat[i][j] = 0.0;
      }
    }
  }
}

void
Dwalk::load(
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
  mat = Matrix(nodesz_+1, Array(nodesz_+1, 0.0)); // **
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
      mat[src_id][dst_id] = w;
      //      if (symm) {
      mat[dst_id][src_id] = w;
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
  lmat  = LabelMatrix(labelsz_+1, LabelArray());
  lcmat = LabelMatrix(labelsz_+1, LabelArray());

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
      lmat[label].push_back(nodesz_);
    }
  }
  ifs.close();
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    unsigned int indx = 0;
    const unsigned int indx_max = lmat[lid].size() - 1;
    for (unsigned int nid = 1; nid <= nodesz_; nid++) {
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
/*
void
Dwalk::cv(
    const unsigned int l_small,
    const unsigned int l_large,
    const std::string& pref_fn,
    const bool map_predict)
{

}
*/

/**
 * apply MAP decision rule and output result
 */
void
Dwalk::decision(
    const Matrix& b,
    const std::string& pref_fn,
    const bool map_predict)
{
  Matrix map(labelsz_+1, Array(nodesz_+1, 0.0));
  Array prob_y(labelsz_+1);
  std::vector<unsigned int> predict(nodesz_+1);

  for (unsigned int i = 1; i <= nodesz_; i++) {
    double sum = 0.0;
    for (unsigned int j = 1; j <= labelsz_; j++) {
      sum += b[j][i];
    }
    map[0][i] = sum;
  }
  for (unsigned int i = 1; i <= nodesz_; i++) {
    for (unsigned int j = 1; j <= labelsz_; j++) {
      map[j][i] = b[j][i] / map[0][i];
      if (isnan(map[j][i])) map[j][i] = 0.0; //**
    }
  }
  // predict
  for (unsigned int i = 1; i <= labelsz_; i++) {
    prob_y[i] = static_cast<double>(lmat[i].size()) / lnodesz_;
  }
  for (unsigned int i = 1; i <= nodesz_; i++) {
    unsigned int predict_tmp = -1;
    double argmax = -1;
    for (unsigned int j = 1; j <= labelsz_; j++) {
      double tmp = map[j][i];
      if (map_predict) tmp *= prob_y[j];
      if (argmax < tmp) {
        predict_tmp = j;
        argmax = tmp;
      }
    }
    predict[i] = predict_tmp;
    assert(1 <= predict[i]);
    assert(predict[i] <= labelsz_);
  }
  // output
  const std::string output_result = pref_fn + ".result";
  const std::string output_weight = pref_fn + ".weight";
  std::ofstream ofs(output_result.c_str());
  ofs.setf(std::ios::fixed, std::ios::floatfield);
  ofs.precision(9);

  for (unsigned int i = 1; i <= nodesz_; i++) {
    ofs << predict[i] << std::endl;
  }
  ofs.close();
  ofs.open(output_weight.c_str());
  for (unsigned int i = 1; i <= nodesz_; i++) {
    for (unsigned int j = 1; j <= labelsz_; j++) {
      if (j != 1) ofs << ' ';
      double tmp = map[j][i];
      if (map_predict) tmp *= prob_y[j];
      ofs << tmp;
    }
    ofs << std::endl;
  }
  ofs.close();
}

/**
 * original Dwalk algorithm:
 *  B_L(q,y) = \frac{
 *    \sum^L \sum \alpha(q) \beta(q)
 *  }{
 *    \sum^L \sum_{q' \in L_y} \alpha(q')
 *  }
 */
void
Dwalk::original(
    const unsigned int bounded_length,
    const std::string& pref_fn,
    const bool denom_flag,
    const bool map_predict)
{
  // calc uset and lset
  NodeSet uset, lset, dummyset;
  std::vector<unsigned int> dummy(nodesz_);
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    lset.insert(lmat[lid].begin(), lmat[lid].end());
  }
  generate(dummy.begin(), dummy.end(), IncrementalGen<unsigned int>(1));
  uset.insert(dummy.begin(), dummy.end());
  set_difference(uset.begin(), uset.end(),
                 lset.begin(), lset.end(),
                 std::inserter(dummyset, dummyset.end()));
  uset.swap(dummyset);
  // calculate forward variable alpha
  // labelsz_, nodesz, mat, lmat
  std::vector<Matrix> alpha_var(bounded_length + 1);
  std::vector<Matrix> beta_var(bounded_length + 1);
  calc_alpha(alpha_var, bounded_length);
  calc_beta (beta_var,  bounded_length);
  // calc bounded dwalks betweeness
  std::cout << "calc bounded dwalks betweeness" << std::endl;
  Matrix b(labelsz_+1, Array(nodesz_+1, 0.0));
  for (unsigned i = 1; i <= labelsz_; i++) {
    if (denom_flag) {
      // calc denom factor
      double denom = 0;
      for (unsigned l = 1; l <= bounded_length; l++) {
        for (unsigned j = 0; j < lmat[i].size(); j++) {
          const unsigned int node_id = lmat[i][j];
          denom += alpha_var[l][i][node_id];
        }
      }
      b[i][0] = 1.0 / denom;
    } else {
      // ignore denom factor
      b[i][0] = 1.0;
    }
    for (unsigned int j = 1; j <= nodesz_; j++) {
      b[i][j] = b[i][0];
    }
  }
  for (unsigned int i = 1; i <= labelsz_; i++) {
    for (unsigned int j = 1; j <= nodesz_; j++) {
      double sum = 0;
      for (unsigned int l = 2; l <= bounded_length; l++) {
        for (unsigned int t = 1; t <= l - 1; t++) {
          sum += alpha_var[t][i][j] * beta_var[l-t][i][j];
          assert(!isnan(sum));
        }
      }
      b[i][j] *= sum;
    }
  }

  decision(b, pref_fn, map_predict);
}

/**
 * fixed algorithm A:
 *  B_L(q,y) = \frac{
 *    \sum^L_l \sum^{l-1}_t \alpha(q) \beta(q)
 *  }{
 *    \sum^L_l \sum^{l-1}_t \alpha(q) \beta(q) + \sum^{L-1}_{t} \alpha(q) \gamma(q)
 *  }
 */
void
Dwalk::fixed_algorithm_A(
    const unsigned int bounded_length,
    const std::string& pref_fn,
    const bool map_predict)
{
  NodeSet uset, lset, dummyset;
  std::vector<unsigned int> dummy(nodesz_);
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    lset.insert(lmat[lid].begin(), lmat[lid].end());
  }
  generate(dummy.begin(), dummy.end(), IncrementalGen<unsigned int>(1));
  uset.insert(dummy.begin(), dummy.end());
  set_difference(uset.begin(), uset.end(),
                 lset.begin(), lset.end(),
                 std::inserter(dummyset, dummyset.end()));
  uset.swap(dummyset);
  std::vector<Matrix> alpha_var(bounded_length + 1);
  std::vector<Matrix> beta_var(bounded_length + 1);
  std::vector<Matrix> gamma_var(bounded_length + 1);
  calc_alpha(alpha_var, bounded_length);
  calc_beta (beta_var,  bounded_length);
  calc_gamma(gamma_var, uset, bounded_length);
  std::cout << "proposed method: " << std::endl;
  std::cout << "calc bounded dwalks betweenness" << std::endl;
  Matrix b(labelsz_+1, Array(nodesz_ + 1, 0.0));
  // calc denom and put b[i][0]
  for (unsigned i = 1; i <= labelsz_; i++) {
    for (unsigned int j = 1; j <= nodesz_; j++) {
      double denom = 0;
      for (unsigned l = 2; l <= bounded_length; l++) {
        for (unsigned t = 1; t <= l - 1; t++) {
          denom += alpha_var[l][i][j] * beta_var[l-t][i][j];
        }
      }
      for (unsigned int t = 1; t <= bounded_length - 1; t++) {
        denom += alpha_var[t][i][j] * gamma_var[bounded_length-t][i][j];
      }
      b[i][j] = 1.0 / denom;
    }
  }
  //
  for (unsigned int i = 1; i <= labelsz_; i++) {
    for (unsigned int j = 1; j <= nodesz_; j++) {
      double sum = 0;
      for (unsigned int l = 2; l <= bounded_length; l++) {
        for (unsigned int t = 1; t <= l - 1; t++) {
          sum += alpha_var[t][i][j] * beta_var[l-t][i][j];
        }
      }
      b[i][j] *= sum;
    }
  }

  decision(b, pref_fn, map_predict);
}

/**
 * fixed algorithm B:
 */
void
Dwalk::fixed_algorithm_B(
    const unsigned int bounded_length,
    const std::string& pref_fn,
    const bool map_predict)
{
  NodeSet uset, lset, dummyset;
  std::vector<unsigned int> dummy(nodesz_);
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    lset.insert(lmat[lid].begin(), lmat[lid].end());
  }
  generate(dummy.begin(), dummy.end(), IncrementalGen<unsigned int>(1));
  uset.insert(dummy.begin(), dummy.end());
  set_difference(uset.begin(), uset.end(),
                  lset.begin(), lset.end(),
                  std::inserter(dummyset, dummyset.end()));
  uset.swap(dummyset);
  std::vector<Matrix> alpha_var(bounded_length + 1);
  std::vector<Matrix> beta_var(bounded_length + 1);
  std::vector<Matrix> gamma_var(bounded_length + 1);
  calc_alpha(alpha_var, bounded_length);
  calc_beta(beta_var,   bounded_length);
  calc_gamma(gamma_var, uset, bounded_length);
  std::cout << "calc bounded dwalks betweenness" << std::endl;
  Matrix b(labelsz_+1, Array(nodesz_ + 1, 0.0));
  // calc denom and put b[i][0]
  for (unsigned int i = 1; i <= labelsz_; i++) {
    for (unsigned int j = 1; j <= nodesz_; j++) {
      double sum = 0;
      for (unsigned int l = 2; l <= bounded_length; l++) {
        for (unsigned int t = 1; t <= l - 1; t++) {
          double tmp = beta_var[l-t][i][j] / (beta_var[l-t][i][j] + gamma_var[l-t][i][j]);
          if (!isnan(tmp)) {
            sum += tmp;
          }
        }
      }
      b[i][j] = sum;
    }
  }
  decision(b, pref_fn, map_predict);
}

} // end of namespace
