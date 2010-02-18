#include "dwalk.h"
namespace gll {

void
Dwalk::showPredict(const std::vector<unsigned int>& predict) const
{
  const unsigned int nodesz = predict.size() - 1;
  for (unsigned int i = 1; i <= nodesz; i++) {
    std::cout << predict[i] << ' ';
  }
  std::cout << std::endl;
}

void
Dwalk::showInfo() const
{
  std::cout << "Number of nodes:           "
            << nodesz_ << std::endl
            << "Number of labels:          "
            << labelsz_ << std::endl
            << "Number of labeled nodes:   "
            << lnodesz_ << std::endl
            << "Number of unlabeled nodes: "
            << unodesz_ << std::endl;
}

const LabelArray
Dwalk::getComplementLabeledNodes (const unsigned int label_id) const
{
  return m_.lcmat_[label_id];
}

const LabelArray
Dwalk::getLabeledNodes (const unsigned int label_id) const
{
  return m_.lmat_[label_id];
}

void
Dwalk::showAlphabeta(
    const unsigned int t,
    const std::string label) const
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(3);
  const unsigned int nodesz = m_.mat_[0].size() - 1;
  std::cout << "forward-backward var:  "
            << label << std::endl
            << "iter:                  "
            << t << std::endl;
  std::cout << "Number of nodes:       "
            << nodesz << std::endl;
  std::cout << "Number of labels:      "
            << labelsz_ << std::endl;
  for (unsigned int j = 1; j <= labelsz_; j++) {
    for (unsigned int i = 1; i <= nodesz; i++) {
      std::cout << m_.mat_[j][i] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Dwalk::showLmat() const
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(3);
  // show lmat
  std::cout << "Number of labels:      "
            << labelsz_ << std::endl;
  for (unsigned int i = 1; i <= labelsz_; i++) {
    std::cout << "label(" << i << "): ";
    for (unsigned int j = 0; j < m_.lmat_[i].size(); j++) {
      std::cout << m_.lmat_[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Dwalk::showBetweenness() const
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(8);
  const unsigned int nodesz = m_.mat_[0].size() - 1;
  for (unsigned int i = 1; i <=labelsz_; i++) {
    for (unsigned int j = 1; j <= nodesz; j++) {
      std::cout << m_.mat_[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Dwalk::showMat() const
{
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(3);
  // show mat
  for (unsigned int i = 1; i <= nodesz_; i++) {
    std::cout << "node " << i << ": ";
    for (unsigned int j = 1; j <= nodesz_; j++) {
      std::cout << m_.mat_[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

}// end of namespace
