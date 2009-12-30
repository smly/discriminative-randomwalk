#AM_CPPFLAGS = $(GTEST_CPPFLAGS)
#AM_CPPFLAGS = $(GTEST_CPPFLAGS) -I$(srcdir)/include
#AM_LDFALGS = $(GTEST_LDFLAGS)

SRC_DIR = ../src
bin_PROGRAMS = test_simple

test_simple_LDADD = $(GTEST_LIBS)
test_simple_CPPFLAGS = -I$(SRC_DIR)
test_simple_SOURCES = test_simple.cc $(SRC_DIR)/dwalk.cc
