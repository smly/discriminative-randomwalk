#AM_CPPFLAGS = $(GTEST_CPPFLAGS)
#AM_CPPFLAGS = $(GTEST_CPPFLAGS) -I$(srcdir)/include
#AM_LDFALGS = $(GTEST_LDFLAGS)

SRC_DIR = ../src
bin_PROGRAMS = test_simple test_realdata

test_simple_LDADD = $(GTEST_LIBS) $(GTEST_LDFLAGS)
test_simple_CPPFLAGS = -I$(SRC_DIR) $(GTEST_CPPFLAGS)
test_simple_SOURCES = test_simple.cc $(SRC_DIR)/dwalk.cc

test_realdata_LDADD = $(GTEST_LIBS) $(GTEST_LDFLAGS)
test_realdata_CPPFLAGS = -I$(SRC_DIR) $(GTEST_CPPFLAGS)
test_realdata_SOURCES = test_realdata.cc $(SRC_DIR)/dwalk.cc