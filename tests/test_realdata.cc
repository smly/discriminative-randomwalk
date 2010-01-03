#include <gtest/gtest.h>
#include <vector>
#include <string>
#include "dwalk.h"

bool exist(const char* filename)
{
  return access(filename,F_OK) == 0;
}

class DwalkTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    graph_fn.push_back(
        std::make_pair("./tests/dataset3/activate-knn9.in",
                       "./tests/dataset3/activate.label"));
    graph_fn.push_back(
        std::make_pair("./tests/dataset3/activate-mat.in",
                       "./tests/dataset3/activate.label"));
  }
  // virtual void TearDown() {}
  std::vector<std::pair<std::string, std::string> > graph_fn;
};

TEST_F(DwalkTest, TestFileExistCheck)
{
  for (unsigned int i = 0; i < graph_fn.size(); i++) {
    const char* graph_filename = graph_fn[i].first.c_str();
    const char* label_filename = graph_fn[i].second.c_str();
    ASSERT_TRUE(exist(graph_filename));
    ASSERT_TRUE(exist(label_filename));
  }
}

TEST_F(DwalkTest, TestFileLoadCheck)
{
  // foreach dataset
  for (unsigned int i = 0; i < graph_fn.size(); i++) {
    const char* graph_filename = graph_fn[i].first.c_str();
    const char* label_filename = graph_fn[i].second.c_str();
    Dwalk d;
    d.load(graph_filename, label_filename, true);// symm
    d.show_info();
    const std::string pref_fn = "output";
    d.go(10, pref_fn, false);
    ASSERT_TRUE(true);
  }
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
