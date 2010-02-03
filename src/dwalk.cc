#include "dwalk.h"

namespace gll {

unsigned int Dwalk::crossValidation(
    const unsigned int num_fold,
    const unsigned int num_from,
    const unsigned int num_to)
{
  return -1;
}

void Dwalk::calcAlpha(
    const unsigned int bounded_length)
{
  Matrix alpha(labelsz_ + 1, Array(nodesz_ + 1, 0.0));
  Matrix alpha_next(labelsz_ + 1, Array(nodesz_ + 1, 0.0));
  // calc \alpha^y(q,1)
  std::cout << "calculate alpha variables: ";
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    for (unsigned int sid = 1; sid <= nodesz_; sid++) {
      const unsigned int nysz = m_.lmat_[lid].size();
      double sum = 0.0;
      for (unsigned int i = 0; i < m_.lmat_[lid].size(); i++) {
        const unsigned int did = m_.lmat_[lid][i];
        sum += (1.0 / nysz) * m_.mat_[did][sid];
        assert(isnan(sum) == 0);
      }
      alpha[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  alpha_[1] = alpha;
  for (unsigned int iter = 2; iter <= bounded_length; iter++) {
    for (unsigned int lid = 1; lid <= labelsz_; lid++) {
      for (unsigned int sid = 1; sid <= nodesz_; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < m_.lcmat_[lid].size(); i++) {
          const unsigned int did = m_.lcmat_[lid][i];
          sum += alpha[lid][did] * m_.mat_[did][sid];
          assert(isnan(sum) == 0);
        }
        alpha_next[lid][sid] = sum;
      }
    }
    alpha = alpha_next;
    alpha_[iter] = alpha;
    std::cout << "." << std::flush;
  }

  std::cout << "done" << std::endl;
}

void Dwalk::calcBeta(
    const unsigned int bounded_length)
{
  Matrix beta(labelsz_+1, Array(nodesz_+1, 0.0));
  Matrix beta_next(labelsz_+1, Array(nodesz_+1, 0.0));
  // calc \beta^y(q,1)
  std::cout << "calculate beta variables:  " << std::flush;
  for (unsigned lid = 1; lid <= labelsz_; lid++) {
    for (unsigned sid = 1; sid <= nodesz_; sid++) {
      double sum = 0.0;
      for (unsigned int i = 0; i < m_.lmat_[lid].size(); i++) {
        const unsigned int did = m_.lmat_[lid][i];
        sum += m_.mat_[sid][did];
        assert(isnan(sum) == 0);
      }
      beta[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  beta_[1] = beta;
  for (unsigned int iter = 2; iter <= bounded_length; iter++) {
    for (unsigned int lid = 1; lid <= labelsz_; lid++) {
      for (unsigned int sid = 1; sid <= nodesz_; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < m_.lcmat_[lid].size(); i++) {
          const unsigned int did = m_.lcmat_[lid][i];
          sum += beta[lid][did] * m_.mat_[sid][did];
          assert(isnan(sum) == 0);
        }
        beta_next[lid][sid] = sum;
      }
    }
    beta = beta_next;
    beta_[iter] = beta;
    std::cout << "." << std::flush;
  }
  std::cout << "done" << std::endl;
}

void Dwalk::calcGamma(
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
        sum += m_.mat_[sid][did];
        assert(isnan(sum) == 0);
      }
      gamma[lid][sid] = sum;
    }
  }
  std::cout << "." << std::flush;
  gamma_[1] = gamma;
  for (unsigned int iter = 2; iter <= bounded_length; iter++) {
    for (unsigned int lid = 1; lid <= labelsz_; lid++) {
      for (unsigned int sid = 1; sid <= nodesz_; sid++) {
        double sum = 0.0;
        for (unsigned int i = 0; i < m_.lcmat_[lid].size(); i++) {
          const unsigned int did = m_.lcmat_[lid][i];
          sum += gamma[lid][did] * m_.mat_[sid][did];
          assert(isnan(sum) == 0);
        }
        gamma_next[lid][sid] = sum;
      }
    }
    gamma = gamma_next;
    gamma_[iter] = gamma;
    std::cout << "." << std::flush;
  }
  std::cout << "done" << std::endl;
}

void
Dwalk::load(
    const char* graph_fn,
    const char* label_fn,
    const bool symm)
{
  m_.load(graph_fn, label_fn, symm);
  nodesz_  = m_.getNodeSize();
  lnodesz_ = m_.getLabeledNodeSize();
  unodesz_ = m_.getUnlabeledNodeSize();
  labelsz_ = m_.getLabelSize();

  assert(m_.lmat_.size() == labelsz_ + 1);
}
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
    prob_y[i] = static_cast<double>(m_.lmat_[i].size()) / lnodesz_;
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

  assert(m_.lmat_.size() == labelsz_ + 1);
  // calc uset and lset
  NodeSet uset, lset, dummyset;
  std::vector<unsigned int> dummy(nodesz_);
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    lset.insert(m_.lmat_[lid].begin(), m_.lmat_[lid].end());
  }
  generate(dummy.begin(), dummy.end(), IncrementalGen<unsigned int>(1));
  uset.insert(dummy.begin(), dummy.end());
  set_difference(uset.begin(), uset.end(),
                 lset.begin(), lset.end(),
                 std::inserter(dummyset, dummyset.end()));
  uset.swap(dummyset);
  // calculate forward variable alpha
  // labelsz_, nodesz, mat, lmat
  alpha_ = std::vector<Matrix>(bounded_length + 1);
  beta_  = std::vector<Matrix>(bounded_length + 1);
  calcAlpha(bounded_length);
  calcBeta (bounded_length);
  // calc bounded dwalks betweeness
  std::cout << "calc bounded dwalks betweeness" << std::endl;
  Matrix b(labelsz_+1, Array(nodesz_+1, 0.0));
  for (unsigned i = 1; i <= labelsz_; i++) {
    if (denom_flag) {
      // calc denom factor
      double denom = 0;
      for (unsigned l = 1; l <= bounded_length; l++) {
        for (unsigned j = 0; j < m_.lmat_[i].size(); j++) {
          const unsigned int node_id = m_.lmat_[i][j];
          denom += alpha_[l][i][node_id];
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
          sum += alpha_[t][i][j] * beta_[l-t][i][j];
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
Dwalk::fixedAlgorithmA(
    const unsigned int bounded_length,
    const std::string& pref_fn,
    const bool map_predict)
{
  NodeSet uset, lset, dummyset;
  std::vector<unsigned int> dummy(nodesz_);
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    lset.insert(m_.lmat_[lid].begin(), m_.lmat_[lid].end());
  }
  generate(dummy.begin(), dummy.end(), IncrementalGen<unsigned int>(1));
  uset.insert(dummy.begin(), dummy.end());
  set_difference(uset.begin(), uset.end(),
                 lset.begin(), lset.end(),
                 std::inserter(dummyset, dummyset.end()));
  uset.swap(dummyset);
  alpha_ = std::vector<Matrix>(bounded_length + 1);
  beta_  = std::vector<Matrix>(bounded_length + 1);
  gamma_ = std::vector<Matrix>(bounded_length + 1);
  calcAlpha(bounded_length);
  calcBeta (bounded_length);
  calcGamma(uset, bounded_length);
  std::cout << "calc bounded dwalks betweenness" << std::endl;
  Matrix b(labelsz_+1, Array(nodesz_ + 1, 0.0));
  // calc denom and put b[i][0]
  for (unsigned i = 1; i <= labelsz_; i++) {
    for (unsigned int j = 1; j <= nodesz_; j++) {
      double denom = 0;
      for (unsigned l = 2; l <= bounded_length; l++) {
        for (unsigned t = 1; t <= l - 1; t++) {
          denom += alpha_[l][i][j] * beta_[l-t][i][j];
        }
      }
      for (unsigned int t = 1; t <= bounded_length - 1; t++) {
        denom += alpha_[t][i][j] * gamma_[bounded_length-t][i][j];
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
          sum += alpha_[t][i][j] * beta_[l-t][i][j];
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
Dwalk::fixedAlgorithmB(
    const unsigned int bounded_length,
    const std::string& pref_fn,
    const bool map_predict)
{
  NodeSet uset, lset, dummyset;
  std::vector<unsigned int> dummy(nodesz_);
  for (unsigned int lid = 1; lid <= labelsz_; lid++) {
    lset.insert(m_.lmat_[lid].begin(), m_.lmat_[lid].end());
  }
  generate(dummy.begin(), dummy.end(), IncrementalGen<unsigned int>(1));
  uset.insert(dummy.begin(), dummy.end());
  set_difference(uset.begin(), uset.end(),
                  lset.begin(), lset.end(),
                  std::inserter(dummyset, dummyset.end()));
  uset.swap(dummyset);
  alpha_ = std::vector<Matrix>(bounded_length + 1);
  beta_  = std::vector<Matrix>(bounded_length + 1);
  gamma_ = std::vector<Matrix>(bounded_length + 1);
  calcAlpha(bounded_length);
  calcBeta(bounded_length);
  calcGamma(uset, bounded_length);
  std::cout << "calc bounded dwalks betweenness" << std::endl;
  Matrix b(labelsz_+1, Array(nodesz_ + 1, 0.0));
  // calc denom and put b[i][0]
  for (unsigned int i = 1; i <= labelsz_; i++) {
    for (unsigned int j = 1; j <= nodesz_; j++) {
      double sum = 0;
      for (unsigned int l = 2; l <= bounded_length; l++) {
        for (unsigned int t = 1; t <= l - 1; t++) {
          double tmp = beta_[l-t][i][j] / (beta_[l-t][i][j] + gamma_[l-t][i][j]);
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
